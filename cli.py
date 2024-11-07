"""tamipami cli module
"""
import argparse
import logging
import tempfile

#import matplotlib.pyplot as plt
#import pandas
#import logomaker

import pam
import degenerate
import fastq
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
                        help='A forward .fastq, .fq, .fastq.gz or .fq.gz file.')
    parser.add_argument('--exp2', '-e2', type=str, default=None,
                        help='A reverse .fastq, .fq, .fastq.gz or .fq.gz file.'),
    parser.add_argument('--outfile', type=str, help='A path to a json file of the results'),
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
