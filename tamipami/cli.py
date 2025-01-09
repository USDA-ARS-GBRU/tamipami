#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""tamipami cli module
"""
import argparse
import logging
import tempfile
import os
import sys
import json

import pandas as pd

from . import pam
from . import degenerate
from . import fastq
from . import tpio
from .config import config


def _logger_setup(logfile: str) -> None:
    """Set up logging to a logfile and the terminal standard out.

    Args:
        logfile: a path to a log file

    Returns:
        None

    """
    try:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%m-%d %H:%M",
            filename=logfile,
            filemode="w",
        )
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter("%(asctime)s: %(levelname)-8s %(message)s")
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger("").addHandler(console)
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
    parser = argparse.ArgumentParser(
        description="TamiPami: a CLI application to parse High throughput PAM/TAM site sequencing data"
    )
    parser.add_argument("--log", help="Log file", default="tamipami.log")
    sub_parsers = parser.add_subparsers(dest="subcommand")

    # create the subparser for the process command
    parser_process = sub_parsers.add_parser(
        "process",
        help='Sub-command "process": used to process FASTQ data into a summarized json output',
    )
    parser_process.add_argument(
        "--cont1",
        "-c",
        type=str,
        required=True,
        help="A forward .fastq, .fq, .fastq.gz or .fq.gz file. .",
    )
    parser_process.add_argument(
        "--cont2",
        "-c2",
        type=str,
        required=True,
        help="A reverse .fastq, .fq, .fastq.gz or .fq.gz file.",
    )
    parser_process.add_argument(
        "--exp1",
        "-e",
        type=str,
        required=True,
        help="A forward .fastq, .fq, .fastq.gz or .fq.gz file.",
    )
    parser_process.add_argument(
        "--exp2",
        "-e2",
        type=str,
        required=True,
        help="A reverse .fastq, .fq, .fastq.gz or .fq.gz file.",
    ),
    parser_process.add_argument(
        "--outfile",
        type=str,
        help="A path to a hdf5 file of the results, if missing binary data will be sent to STDOUT.",
    ),
    parser_process.add_argument(
        "--library",
        type=str,
        choices=["RTW554", "RTW555", "RTW572", "RTW574"],
        required=False,
        help="The Addgene library pool. For custom pools use the --spacer and --orientation flags",
    )
    parser_process.add_argument(
        "--spacer",
        type=str,
        required=False,
        help="The spacer sequence for the guide RNA. Not needed if ---library is used",
    )
    parser_process.add_argument(
        "--orientation",
        type=str,
        choices=["3prime", "5prime"],
        required=False,
        help="the side of the spacer the PAM/TAM is on",
    )
    parser_process.add_argument(
        "--length",
        choices=range(3, 9),
        metavar="[3-8]",
        default=6,
        type=int,
        help=" The length of the PAM or TAM sequences",
    )

    # create a subparser for the predict command to provide the degenerate sequences and output data  given a cutoff and length
    parser_predict = sub_parsers.add_parser(
        "predict",
        help='Subcommand "predict" use to predict PAMs/TAMs and summary data for a selected length and cutoff value',
    )
    parser_predict.add_argument(
        "--input",
        type=str,
        required=False,
        help="An file containing the data from a TamiPami process run or downloaded data from the web ap, if not input is provided STDIN will be assumed",
    )
    parser_predict.add_argument(
        "--cutoff",
        type=str,
        required=False,
        help="""A json string containing the kmer lengths and the Zscore cutoff values above which kmers are considered part of the PAM/TAM.
                               If no cutoff is provided it will be automatically calculated using univariate k means clustering. example input : '{'4': 0.7, '5': 1.35}'
                               """,
    )
    parser_predict.add_argument(
        "--predict_out",
        type=str,
        required=True,
        help="A file directory containing the PAM/TAM and degenerate sequences identified",
    )

    return parser


spacer_dict = config["spacer_dict"]


def parse_lib(args: argparse.Namespace) -> tuple[str]:
    if args.library:
        spacer = config["spacer_dict"][args.library]["spacer"]
        orientation = config["spacer_dict"][args.library]["orientation"]
    else:
        spacer = args.spacer
        orientation = args.orientation
    return spacer, orientation


def process(args: argparse.Namespace = None) -> tuple[pam.pamSeqExp, dict]:
    try:
        with tempfile.TemporaryDirectory() as datadir:
            spacer, orientation = parse_lib(args)
            logging.info("spacer: {}".format(spacer))
            logging.info("orientation: {}".format(orientation))
            run_summ = {"cont": {}, "exp": {}}
            cont_raw, run_summ["cont"]["tot"], run_summ["cont"]["targets"] = (
                fastq.process(
                    fastq=args.cont1,
                    fastq2=args.cont2,
                    pamlen=args.length,
                    mergedfile=os.path.join(datadir, "cont_merged.fastq.gz"),
                    spacer=spacer,
                    orientation=orientation,
                )
            )
            exp_raw, run_summ["exp"]["tot"], run_summ["exp"]["targets"] = fastq.process(
                fastq=args.exp1,
                fastq2=args.exp2,
                pamlen=args.length,
                mergedfile=os.path.join(datadir, "exp_merged.fastq.gz"),
                spacer=spacer,
                orientation=orientation,
            )
            pamexpobj = pam.pamSeqExp(ctl=cont_raw, exp=exp_raw, position=orientation)
            return pamexpobj, run_summ
    except Exception as e:
        logging.error("Error processing the fastq files")
        raise e


import os
import logging


def _create_directories(outdir: str, subdir_list: list) -> str:
    """Create output dir
    Args:
        outdir: the path to the outputdir
        subdir_list: a list of subdirs to create

    Returns:
    str: the mane of the outdir filepath
    """

    try:
        # Check if the main directory exists
        if os.path.exists(outdir):
            logging.error(f"Directory '{outdir}' already exists.")
            raise FileExistsError(f"Directory '{outdir}' already exists.")

        # Create the main directory
        os.makedirs(outdir)

        # Create the subdirectories
        for subdir in subdir_list:
            subdir_path = os.path.join(outdir, subdir)
            os.makedirs(subdir_path, exist_ok=True)

    except FileExistsError as fee:
        logging.error(
            f"The output directory already exists please remove it or specify a new location {fee}"
        )
        raise

    except OSError as e:
        logging.error(f"Error creating directory '{outdir}': {e}")
        raise

    return outdir


def _mk_outputs(pamseqobj: pam.pamSeqExp, outdir: str, lenval: int, cutoff: float):
    try:
        df = pamseqobj.multikmerdict[lenval]
        pamseqobj.make_logo(
            length=lenval,
            cutoff=cutoff,
            type="zscore",
            above=True,
            filename=os.path.join(
                outdir, str(lenval), "logo." + config["logo_file_type"]
            ),
        )
        df.to_csv(os.path.join(outdir, str(lenval), "data.csv"), index=False)
        tpio.histogram_plot(
            df,
            maxbins=config["histogram_bins"],
            cutoff=cutoff,
            filename=os.path.join(outdir, str(lenval), "histogram.html"),
        )
        filtered_df = df[df["zscore"] >= cutoff]
        dseqs = degenerate.seqs_to_degenerates(filtered_df["kmers"].tolist())
        degen_df = pd.DataFrame(dseqs, columns=["PAM/TAM site"])
        degen_df.to_csv(
            os.path.join(outdir, str(lenval), "degenerate_pam_tam.csv"), index=False
        )
    except Exception as e:
        logging.error(
            "could not generate summarized export data for kmern length {}".format(
                lenval
            )
        )
        raise e


def predict(args: argparse.Namespace = None) -> None:
    if args.input:
        data = tpio.import_hdf(from_buffer=False, filename=args.input)
    else:
        data = tpio.import_hdf(from_buffer=True)
    pamseqobj = pam.pamSeqExp(multikmerdict=data)
    kmdkeys = list(pamseqobj.multikmerdict.keys())
    print(kmdkeys)
    if args.cutoff:
        cutoff = json.loads(args.cutoff)
        assert set(cutoff.keys()) == set(
            kmdkeys
        ), "a length listed in the --cutoff input does not match the lengths available in "
    else:
        cutoff = {}
        for lenval in kmdkeys:
            cutoff[lenval] = pamseqobj.find_breakpoint(length=lenval)
    pd_out = _create_directories(
        args.predict_out, [str(i) for i in list(pamseqobj.multikmerdict.keys())]
    )
    for lenval in kmdkeys:
        _mk_outputs(
            pamseqobj=pamseqobj, outdir=pd_out, lenval=lenval, cutoff=cutoff[lenval]
        )


def main(args: argparse.Namespace = None) -> None:
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
    if args.subcommand == "process":
        pamexpobj, run_summ = process(args)
        logging.info(
            "Control data: {} reads merged, {} contained targets.".format(
                run_summ["cont"]["tot"], run_summ["cont"]["targets"]
            )
        )
        logging.info(
            "Experimental data: {} reads merged, {} contained targets.".format(
                run_summ["exp"]["tot"], run_summ["exp"]["targets"]
            )
        )
        logging.info("Writing results to {}".format(args.outfile))
        if args.outfile:
            tpio.export_hdf(
                multikmerdict=pamexpobj.multikmerdict,
                to_buffer=False,
                filename=args.outfile,
            )
        else:
            sys.stdout.buffer.write(
                tpio.export_hdf(
                    multikmerdict=pamexpobj.multikmerdict,
                    to_buffer=True,
                    filename="tamipami.h5",
                )
            )
    elif args.subcommand == "predict":
        predict(args)


if __name__ == "__main__":
    main()
