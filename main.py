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
from collections import Counter
import statistics
from typing import Dict, List, Tuple, Literal

import matplotlib.pyplot as plt
import pandas
from Bio import SeqIO
import scipy.stats
from skbio.stats import composition
from statsmodels.stats.multitest import fdrcorrection
import logomaker


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
    parser = argparse.ArgumentParser(description='HT-TAMDA: a script to parse High throughput PAM sequences')
    parser.add_argument('--cont1', '-c', type=str, required=False,
                        help='A forward .fastq, .fq, .fastq.gz or .fq.gz file. .')
    parser.add_argument('--cont2', '-c2', type=str, required=False,
                        help='A reverse .fastq, .fq, .fastq.gz or .fq.gz file.')
    parser.add_argument('--exp1', '-e', type=str, required=False,
                        help='A forward .fastq, .fq, .fastq.gz or .fq.gz file. .')
    parser.add_argument('--exp2', '-e2', type=str, default=None,
                        help='A reverse .fastq, .fq, .fastq.gz or .fq.gz file.')
    parser.add_argument('--library', type=str, choices=["RTW544", "RTW555", "RTW572", "RTW574"],
                        help='The Addgene library pool. For custom pools use the --spacer and --orientation flags'
                        )
    parser.add_argument('--spacer', type=str,  help='The spacee sequence for the guide RNA. Not needed if ---library is used' )
    parser.add_argument('--orientation', type=str, choices=["3prime", "5prime"], help='the side of the sacer the PAM/TAM is on')
    parser.add_argument('--log', help="Log file", default="HT-TAMDA.log")
    parser.add_argument('--length', choices=range(1, 11), metavar="[1-10]",
                        help=" The length of the PAM or TAM sequences", default=4, type=int)
    # parser.add_argument('--threads', help="Number of processor threads to use.", type=int, default=1)
    return parser


spacer_dict = {'RTW572': 'GGAATCCCTTCTGCAGCACCTGG',
               'RTW554': 'GGGCACGGGCAGCTTGCCGG',
               'RTW555': 'GTCGCCCTCGAACTTCACCT',
               'RTW574': 'CTGATGGTCCATGTCTGTTACTC'}

# spacer_dict = {'RTW572': {'spacer': 'GGAATCCCTTCTGCAGCACCTGG','orientation': '3prime'},
#               'RTW554': {'spacer': 'GGGCACGGGCAGCTTGCCGG', 'orientation': '3prime'},
#               'RTW555': {'spacer': 'GTCGCCCTCGAACTTCACCT', 'orientation': '5prime'},
#               'RTW574': {'spacer': 'CTGATGGTCCATGTCTGTTACTC', 'orientation': '5prime'}


def iterate_kmer(k: int) -> Dict:
    """Create a list of kmers.

    Args:
        k : The kmer  with a size representing the length of the PAM or TAM sequence, typically 4 or 5.

    Returns:
        A dict with lexographically sorted kmers as keys and 0 as values.
    """

    bases = ['A', 'C', 'T', 'G']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    return dict.fromkeys(kmers, 0)


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
                      'out=' + outfile,
                      'ecco=t',
                      'mix=f']
        p3 = subprocess.run(parameters, stderr=subprocess.PIPE)
        return p3.stderr.decode('utf-8')
    except RuntimeError:
        logging.warning("could not perform read merging with bbmerge")


def count_pam(spacer: str, fastq: str, pamlen: int, orientation: str) -> Dict:
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
    for record in SeqIO.parse(fastq, "fastq"):
        seqstr = str(record.seq)
        result = re.search(spacer, seqstr)
        if result:
            if orientation == '5prime':
                spacerstart = result.start()
                pamstart = spacerstart - int(pamlen)
                pamseq = seqstr[pamstart:spacerstart]
            if orientation== '3prime':
                pass
            if pamseq in kmer_dict:
                kmer_dict[pamseq] += 1
    return kmer_dict


def process(fastq: str, fastq2: str, pamlen: int, tempdir: str, spacer: str, orientation: str) -> Tuple[Dict, List]:
    """A function to merge the reads and count the TAM/PAM sequences
    
    Args:
        fastq (str): 
        fastq2 (str): 
        pamlen (int: 
        tempdir: 
        spacer: 

    Returns:
        A Dictionay contining counts of every  PAM/TAM in a sample
    """
    mergedfile = os.path.join(tempdir, "merged.fastq")
    logging.info("Merging reads.")
    stdout = merge_reads(fastq=fastq, fastq2=fastq2, outfile=mergedfile)
    logging.info(stdout)
    logging.info("Counting PAM/TAM sites.")
    refcount = count_pam(pamlen=pamlen, spacer=spacer, fastq=mergedfile, orientation=orientation)
    logging.info("Runing multiplicitive replacment in case there are zeros")
    refmc = composition.multiplicative_replacement(list(refcount.values()))
    logging.info("Calculating a center log ratio to change from Aitchison geometry to the real number space")
    refclr = composition.clr(refmc)
    ref_n = check_N(list(refcount.values()))
    logging.info("Poisson noise is expected to be {:5.1f} % of total".format(ref_n * 100))
    return refcount, refclr


def check_N(vect: List) -> float:
    """ Estimate the fraction of SD expected to come from shot noise

    The fraction of the SD expected to be attributable to Poisson or shot noise.
      A large value indicates that more reads or a shorter length value are required.
      
    Args:
        vect (list): A list of  kmer counts
        
    Returns:
        float: The fraction of SD attributable Poisson Noise
    """
    return math.sqrt(statistics.mean(vect)) / statistics.stdev(vect)


def make_df(cont_raw: Dict, cont_clr: List, exp_raw: Dict, exp_clr: List) -> pandas.DataFrame:
    """ Make a data frame and estimate Z scores and BH adjusted p-values for the clr(control) -clr(experimental
    Args:
        cont_raw: A Dictionary of PAM/TAM counts for the control library
        cont_clr: A List of  center Log Ratio Transformed Values for the control library
        exp_raw: A Dictionary of PAM/TAM counts for the experimental library
        exp_clr: A List of  center Log Ratio transformed Values for the experimental library

    Returns:
        A Pandas dataframe with rows sorted by their deviation from controls
    """
    data = {"seqs": list(cont_raw.keys()),
            "cont_raw": list(cont_raw.values()),
            "exp_raw": list(exp_raw.values()),
            "cont_clr": list(cont_clr),
            "exp_clr": list(exp_clr)}
    df = pandas.DataFrame.from_dict(data)
    df["diff"] = df["cont_clr"] - df["exp_clr"]
    df["zscore"] = (df['diff'] - df['diff'].mean())/ df['diff'].std()
    # ERROR in use of CDF!!!!
    df["pvalue"] = scipy.stats.norm.pdf(df['zscore'])
    df["signifigant"], df["p_adjust_BH"] = fdrcorrection(df["pvalue"], alpha=0.05)
    print(df)
    df = df.sort_values(by=["p_adjust_BH"])
    return df

def make_logo(df: pandas.DataFrame, padjust: float, filename: str, ) -> None:
    """
    Takes a pandas dataframe and saves a Sequence motif logo of the signifigant

    Args:
        df: A data frame genrated by the make_df function
        padjust: the threshold for adjusted pvalues to include in the motif
        filename:  a path for the pdf, jpg or png of figure
    Returns:
        None
    """
    try:
        df_filtered = df[df["p_adjust_BH"] <= padjust]
        print("df_filtered")
        print(df_filtered)
        prob_df = logomaker.alignment_to_matrix(sequences=df_filtered['seqs'], to_type = 'probability', pseudocount = 0)
        logo_fig = logomaker.Logo(prob_df, color_scheme='classic')
        plt.savefig(filename)
    except Exception as e:
        logging.warning( "could not generate sequence motif graphic")
        raise e


def main(args=None):
    """Run the TAM/PAM identif≈ication workflow

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
    logging.info("Processing control reads")
    tempdir = tempfile.mkdtemp()
    cont_raw, cont_clr = process(fastq=args.cont1, fastq2=args.cont2, pamlen=args.length, tempdir=tempdir,
                               spacer=spacer_dict[args.library], orientation='5prime')
    tempdir2 = tempfile.mkdtemp()
    logging.info("Processing experimental reads")
    exp_raw, exp_clr = process(fastq=args.exp1, fastq2=args.exp2, pamlen=args.length, tempdir=tempdir2,
                               spacer=spacer_dict[args.library], orientation='5prime')
    df = make_df(cont_raw=cont_raw, cont_clr=cont_clr, exp_raw=exp_raw, exp_clr=exp_clr)
    make_logo(df=df, padjust=0.05, filename="logo.pdf" )


if __name__ == "__main__":
    main()
