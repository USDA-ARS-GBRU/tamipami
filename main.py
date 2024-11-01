import itertools
# -*- coding: utf-8 -*-
# #!/usr/bin/local python
"""
A script to process High throghput PAM/TAM sequencing libraries

"""
import math
import subprocess
import argparse
import logging
import tempfile
import os
import re
import itertools
import statistics
import gzip

#import matplotlib.pyplot as plt
#import pandas
from Bio import SeqIO
#import logomaker

import pam
import degenerate

from config import config


def _logger_setup(logfile: str) -> None:
    """Set up logging to a logfile and the terminal standard out.

    Args:
        logfile: a path to a log file

    Returns:
        None

    """
    try:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename=logfile,
                            filemode='w')
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger('').addHandler(console)
    except Exception as e:
        print("An error occurred setting up logging")
        raise e


def myparser() -> argparse.ArgumentParser:
    """Set up logging to a logfile and the terminal standard out.
    Args:
        None
    Returns:
        A filled ArgumentParser object
    """
    parser = argparse.ArgumentParser(description='tamipami: a script to parse High throughput PAM sequences')
    parser.add_argument('--cont1', '-c', type=str, required=False,
                        help='A forward .fastq, .fq, .fastq.gz or .fq.gz file. .')
    parser.add_argument('--cont2', '-c2', type=str, required=False,
                        help='A reverse .fastq, .fq, .fastq.gz or .fq.gz file.')
    parser.add_argument('--exp1', '-e', type=str, required=False,
                        help='A forward .fastq, .fq, .fastq.gz or .fq.gz file. .')
    parser.add_argument('--exp2', '-e2', type=str, default=None,
                        help='A reverse .fastq, .fq, .fastq.gz or .fq.gz file.')
    parser.add_argument('--library', type=str, choices=["RTW554", "RTW555", "RTW572", "RTW574"],
                        help='The Addgene library pool. For custom pools use the --spacer and --orientation flags'
                        )
    parser.add_argument('--spacer', type=str,  help='The spacee sequence for the guide RNA. Not needed if ---library is used' )
    parser.add_argument('--orientation', type=str, choices=["3prime", "5prime"], help='the side of the sacer the PAM/TAM is on')
    parser.add_argument('--log', help="Log file", default="tamipami.log")
    parser.add_argument('--length', choices=range(1, 11), metavar="[1-10]",
                        help=" The length of the PAM or TAM sequences", default=4, type=int)
    # parser.add_argument('--threads', help="Number of processor threads to use.", type=int, default=1)
    return parser

spacer_dict  = config["spacer_dict"]


def iterate_kmer(k: int) -> dict[str,int]:
    """
    Create a dictionary of all possible kmers of length k. 
    Args:
        k: The length of the kmer.      
    Returns:
        A dictionary with lexographically sorted kmers as keys and 0 as values.
    """
    bases = ['A', 'C', 'T', 'G']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    return {kmer: 0 for kmer in sorted(kmers)}

def merge_reads(fastq, fastq2, outfile) -> str:
    """ Merge Reads and return error-corrected merged fastq

    Args:
        fastq: The path to the forward fastq file.
        fastq2: The path to reverse fastq file.
        outfile: The path to the merged fastq output.

    Returns:
        The stdout from BBmerge
    """

    try:
        parameters = ['bbmerge.sh', 'in=' + fastq,
                      'in2=' + fastq2,
                      'out=' + outfile,]
        parameters.extend(config["bbmerge"])
        p3 = subprocess.run(parameters, stderr=subprocess.PIPE)
        return p3.stderr.decode('utf-8')
    except RuntimeError:
        logging.warning("could not perform read merging with bbmerge")


def count_pam(spacer: str, fastq: str, pamlen: int, orientation: str) -> dict[str, int]:
    """Create dictionary of PAM/TAM counts

    Args:
        spacer: The DNA sequence matching the guide sequence
        fastq: The  path to a merged fastq file
        pamlen: The legnth of the PAM/TAM sequences
        orientation: 5prime (5'-PAM-Spacer-3') or 3prime (5'-Spacer-PAM-3')

    Returns:
        A dictionary of kmer counts
    """
    kmer_dict = iterate_kmer(pamlen)
    with gzip.open(fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seqstr = str(record.seq)
            result = re.search(spacer, seqstr)
            if result:
                if orientation == '5prime':
                    spacerstart = result.start()
                    pamstart = spacerstart - int(pamlen)
                    pamseq = seqstr[pamstart:spacerstart]
                elif orientation== '3prime':
                    spacerend = result.end()
                    pamend = spacerend + int(pamlen)
                    pamseq = seqstr[spacerend:pamend]
                if pamseq in kmer_dict:
                    kmer_dict[pamseq] += 1
    return kmer_dict


def process(fastq: str, fastq2: str, pamlen: int, tempdir: str, spacer: str, orientation: str) -> tuple[dict, list]:
    """A function to merge the reads and count the TAM/PAM sequences

    Args:
        fastq (str): 
        fastq2 (str): 
        pamlen (int: 
        tempdir: 
        spacer: 

    Returns:
    A  Dictionay contining counts of every  PAM/TAM in a sample
"""
    mergedfile = os.path.join(tempdir, "merged.fastq.gz")
    logging.info("Merging reads.")
    stdout = merge_reads(fastq=fastq, fastq2=fastq2, outfile=mergedfile)
    logging.info(stdout)
    logging.info("Counting PAM/TAM sites.")
    pamcount = count_pam(pamlen=pamlen, spacer=spacer, fastq=mergedfile, orientation=orientation)
    ref_n = check_N(list(pamcount.values()))
    logging.info("Poisson noise is expected to be {:5.1f} % of total at N={}".format(ref_n * 100, pamlen))
    return pamcount


def check_N(vect: list) -> float:
    """ Estimate the fraction of SD expected to come from shot noise
        The fraction of the SD expected to be attributable to Poisson or shot noise.
        A large value indicates that more reads or a shorter length value are required.  
    Args:
        vect (list): A list of kmer counts  
    Returns:
        float: The fraction of SD attributable Poisson Noise
    """
    if not vect:
        return 0.0
    return math.sqrt(statistics.mean(vect)) / statistics.stdev(vect)

# def make_logo(df: pandas.DataFrame, padjust: float, filename: str, ) -> None:
#     """
#     Takes a pandas dataframe and saves a Sequence motif logo of the signifigant

#     Args:
#         df: A data frame genrated by the make_df function
#         padjust: the threshold for adjusted pvalues to include in the motif
#         filename:  a path for the pdf, jpg or png of figure
#     Returns:
#         None
#     """
#     try:
#         df_filtered = df[df["p_adjust_BH"] <= padjust]
#         prob_df = logomaker.alignment_to_matrix(sequences=df_filtered['seqs'], to_type = 'probability', pseudocount = 0)
#         logo_fig = logomaker.Logo(prob_df, color_scheme='classic')
#         plt.savefig(filename)
#     except Exception as e:
#         logging.warning( "could not generate sequence motif graphic")
#         raise e


def main(args: argparse.Namespace=None) -> None:
    """Run the TAM/PAM identification workflow
    Args:
        args: arg parse args to pass in (used for unit testing)
    Returns:
        None
    """
    parser = myparser()
    if not args:
        args = parser.parse_args()
    _logger_setup(args.log)
    logging.info("Begin processing PAM/TAM sequencing libraries")
    logging.info(args)
    if args.library:
        spacer = spacer_dict[args.library]["spacer"]
        orientation = spacer_dict[args.library]["orientation"]
    else:
        spacer = args.spacer
        orientation = args.orientation
    logging.info("Processing control reads")
    cont_raw = process(fastq=args.cont1, fastq2=args.cont2, pamlen=args.length, tempdir=tempfile.mkdtemp(),
                               spacer=spacer, orientation=orientation)
    logging.info("Processing experimental reads")
    exp_raw = process(fastq=args.exp1, fastq2=args.exp2, pamlen=args.length, tempdir=tempfile.mkdtemp(),
                               spacer=spacer, orientation=orientation)
    #df = make_df(cont_raw=cont_raw, cont_clr=cont_clr, exp_raw=exp_raw, exp_clr=exp_clr)
    logging.info("creating a PamSeqExp object")
    pamexpobj = pam.pamSeqExp(ctl=cont_raw, exp=exp_raw,position=orientation)
    pamexpobj.plot_kmer_summary(attribute='zscore', save_path="summaryplot.pdf")
    breakpoint = pamexpobj.find_breakpoint(length=args.length, type='diff')
    largestk = pamexpobj.multikmerdict[int(args.length)]
    cleaved_seqs = largestk[largestk['diff'] < breakpoint]
    degenerate_seqs = degenerate.seqs_to_degenerates(list(cleaved_seqs['kmers']))
    logging.info("the Pam sequences are {}".format(degenerate_seqs))


if __name__ == "__main__":
    main()
