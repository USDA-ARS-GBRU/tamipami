#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""cli: a TamiPami command line interfaces application module"""
import argparse
import logging
import tempfile
import os
import sys
import multiprocessing
from pathlib import Path

import pandas as pd

from tamipami import pam
from tamipami import degenerate
from tamipami import fastq
from tamipami import tpio
from tamipami.cli_utils import cutoff_arg_validator
from tamipami.config import config
from tamipami._version import __version__


def _logger_setup(logfile: str) -> None:
    """Set up logging configuration to output logs to a specified file and the console.

    This function configures the logging system to write logs to a specified
    logfile and also outputs logs of level INFO or higher to the console. The
    log level, format, and date format can be customized using environment
    variables.

    Args:
        logfile: The path to the log file where logs will be written.

    Raises:
        FileNotFoundError: If the directory for the logfile does not exist.
    """
    logdir = os.path.dirname(logfile)
    if logdir and not os.path.exists(logdir):
        os.makedirs(logdir, exist_ok=True)

    log_level = os.getenv("LOG_LEVEL", "DEBUG")
    log_format = os.getenv(
        "LOG_FORMAT", "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    )
    date_format = os.getenv("DATE_FORMAT", "%m-%d %H:%M")

    logging.basicConfig(
        level=getattr(logging, log_level),
        format=log_format,
        datefmt=date_format,
        filename=logfile,
        filemode="w",
    )
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console: logging.StreamHandler = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter: logging.Formatter = logging.Formatter(
        "%(asctime)s: %(levelname)-8s %(message)s"
    )
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger("").addHandler(console)


def myparser() -> argparse.ArgumentParser:
    """
    Creates and configures an ArgumentParser for the TamiPami CLI application.

    This function sets up the argument parser with subcommands and options
    for processing and predicting PAM/TAM sequences from high throughput
    sequencing data. It includes options for specifying input files, output
    locations, and various parameters related to the sequencing data and
    analysis.

    Returns:
        argparse.ArgumentParser: A configured ArgumentParser object with
        subcommands and arguments for the TamiPami CLI.
    """
    parser = argparse.ArgumentParser(
        description="TamiPami: a CLI application to parse High throughput PAM/TAM site sequencing data"
    )
    parser.add_argument("--log", help="Log file", default="tamipami.log")
    sub_parsers = parser.add_subparsers(dest="subcommand")

    parser.add_argument("-V", "--version", action="version", version=__version__)

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
        required=False,
        help="Optional reverse .fastq, .fq, .fastq.gz or .fq.gz file for paired-end runs. Provide both --cont2 and --exp2, or omit both for single-end.",
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
        required=False,
        help="Optional reverse .fastq, .fq, .fastq.gz or .fq.gz file for paired-end runs. Provide both --cont2 and --exp2, or omit both for single-end.",
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
        type=cutoff_arg_validator,
        required=True,
        help=(
            "Cutoff thresholds as a JSON dictionary with integer keys [3-8] and numeric values. "
            "Example: '{\"3\": 0.7, \"4\": 0.85, \"5\": 0.93}'. "
            "Keys must be integers between 3 and 8. Values must be numbers."
        ),
    )
    parser_predict.add_argument(
        "--predict_out",
        type=str,
        required=True,
        help="A file directory containing the PAM/TAM and degenerate sequences identified",
    )

    return parser


spacer_dict = config["spacer_dict"]


def parse_lib(args: argparse.Namespace) -> tuple[str, str]:
    """
    Parse library-specific configuration for spacer and orientation.

    Args:
        args (argparse.Namespace): The command-line arguments containing library,
        spacer, and orientation information.

    Returns:
        tuple[str, str]: A tuple containing the spacer and orientation values.
        If a library is specified, these values are retrieved from the configuration;
        otherwise, they default to the provided arguments.
    """
    if args.library:
        spacer = config["spacer_dict"].get(args.library, {}).get("spacer", args.spacer)
        orientation = (
            config["spacer_dict"]
            .get(args.library, {})
            .get("orientation", args.orientation)
        )
    else:
        spacer = args.spacer
        orientation = args.orientation
    return spacer, orientation

def process_fastq_wrapper(args_tuple):
    """
    Wrapper function to process FASTQ files using provided arguments with multiprocessing.

    Args:
        args_tuple (tuple): Contains (fastq1, fastq2, spacer, orientation, pamlen).

    Returns:
        tuple: Result from fastq.process, including kmer counts, total reads, and guide detections.
    """
    fastq1, fastq2, spacer, orientation, pamlen = args_tuple
    return fastq.process(
        fastq=fastq1,
        fastq2=fastq2,
        pamlen=pamlen,
        spacer=spacer,
        orientation=orientation,
    )

def process(args: argparse.Namespace = None) -> tuple[pam.pamSeqExp, dict]:
    """
    Process control and experimental FASTQ files to generate a pamSeqExp object and run summary.

    Args:
        args (argparse.Namespace, optional): Command-line arguments specifying input files, spacer, orientation, and kmer length.

    Returns:
        tuple[pam.pamSeqExp, dict]: A tuple containing the pamSeqExp object and a run summary dictionary.

    Raises:
        FileNotFoundError: If an input file is not found.
        ValueError: If input arguments are invalid.
        Exception: For unexpected errors during processing.
    """
    try:
        with tempfile.TemporaryDirectory() as datadir:
            os.chmod(datadir, 0o700)
            spacer, orientation = parse_lib(args)
            logging.info(f"spacer: {spacer}")
            logging.info(f"orientation: {orientation}")
            run_summ = {"cont": {}, "exp": {}}

            # Validate paired-end vs single-end inputs
            if bool(args.cont2) != bool(args.exp2):
                raise ValueError(
                    "Invalid input: provide both --cont2 and --exp2 for paired-end processing, or omit both for single-end."
                )
            paired_end = bool(args.cont2) and bool(args.exp2)
            if paired_end:
                logging.info("Paired-end mode enabled.")
            else:
                logging.info("Single-end mode enabled.")

            # Prepare arguments for both calls
            tasks = [
                (args.cont1, args.cont2, spacer, orientation, args.length),
                (args.exp1, args.exp2, spacer, orientation, args.length),
            ]

            with multiprocessing.Pool(processes=2) as pool:
                results = pool.map(process_fastq_wrapper, tasks)

            (cont_raw, run_summ["cont"]["tot"], run_summ["cont"]["targets"]) = results[0]
            (exp_raw, run_summ["exp"]["tot"], run_summ["exp"]["targets"]) = results[1]

            pamexpobj = pam.pamSeqExp(ctl=cont_raw, exp=exp_raw, position=orientation)
            return pamexpobj, run_summ
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        raise
    except ValueError as e:
        logging.error(f"Value error: {e}")
        raise
    except Exception as e:
        logging.error("Unexpected error processing the fastq files")
        raise


import os
import logging


def _create_directories(outdir: str, subdir_list: list) -> str:
    """Create a main output directory and specified subdirectories.

    Args:
        outdir (str): The path to the main output directory.
        subdir_list (list): A list of subdirectories to create within the main directory.

    Raises:
        FileExistsError: If the output directory already exists.
        OSError: If there is an error creating the directory.

    Logs a warning if the subdir_list is empty and logs errors for exceptions.
    """
    try:
        # Create the main directory
        Path(outdir).mkdir(parents=True, exist_ok=True)

        # Check if subdir_list is empty
        if not subdir_list:
            logging.warning(
                "The subdir_list is empty. No subdirectories will be created."
            )

        # Create the subdirectories
        for subdir in subdir_list:
            subdir_path = Path(outdir) / subdir
            subdir_path.mkdir(parents=True, exist_ok=True)

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
    """
    Generates and saves outputs for a given PAM/TAM sequence experiment.

    This function processes kmer data from a pamSeqExp object, generating a sequence
    motif logo, histogram plot, and degenerate sequence file for a specified kmer length
    and cutoff. It saves these outputs to a specified directory.

    Args:
        pamseqobj (pam.pamSeqExp): The pamSeqExp object containing multikmer data.
        outdir (str): The directory path where outputs will be saved.
        lenval (int): The kmer length to process.
        cutoff (float): The z-score cutoff for filtering kmers.

    Raises:
        ValueError: If `outdir` is not a valid directory, `lenval` is not a positive
        integer.
        KeyError: If the specified kmer length is not found in multikmerdict.
        FileNotFoundError: If a file cannot be saved in the specified directory.
        Exception: For any unexpected errors during processing.
    """
    if not isinstance(outdir, str) or not os.path.isdir(outdir):
        raise ValueError("`outdir` must be a valid directory path.")
    if not isinstance(lenval, int) or lenval <= 0:
        raise ValueError("`lenval` must be a positive integer.")
    try:
        df = pamseqobj.multikmerdict[lenval]
        pamseqobj.make_logo(
            length=lenval,
            cutoff=cutoff,
            above=True,
            filename=os.path.join(
                outdir, str(lenval), "logo." + config["logo_file_type"]
            ),
        )
        os.makedirs(os.path.join(outdir, str(lenval)), exist_ok=True)
        df.to_csv(os.path.join(outdir, str(lenval), "data.csv"), index=False)
        maxbins = config.get("histogram_bins", 10)  # Default to 10 if not found
        tpio.histogram_plot(
            df,
            maxbins=maxbins,
            cutoff=cutoff,
            filename=os.path.join(outdir, str(lenval), "histogram.html"),
        )
        filtered_df = df[df["zscore"] >= cutoff]
        dseqs = degenerate.seqs_to_degenerates(filtered_df["kmers"].tolist())
        degen_df = pd.DataFrame(dseqs, columns=["PAM/TAM site"])
        degen_df.to_csv(
            os.path.join(outdir, str(lenval), "degenerate_pam_tam.csv"), index=False
        )
    except KeyError as e:
        logging.error(
            f"KeyError: Could not find kmer length {lenval} in multikmerdict."
        )
        raise e
    except FileNotFoundError as e:
        logging.error(f"FileNotFoundError: Could not save file in directory {outdir}.")
        raise e
    except Exception as e:
        logging.error(f"Unexpected error occurred for kmer length {lenval}: {str(e)}")
        raise e


def predict(args: argparse.Namespace = None) -> None:
    """
    Predicts TAM/PAM sites using input data from an HDF file or stream.

    This function imports multikmer data, applies a cutoff value to filter
    the data, and generates output directories and files for each kmer length.
    If no cutoff is provided, it calculates cutoff breakpoints for each kmer length.

    Args:
        args (argparse.Namespace): Command-line arguments containing input
        file paths and cutoff values.

    Raises:
        AssertionError: If the provided cutoff lengths do not match the
        available kmer lengths in the data.
    """
    data = tpio.import_hdf(from_buffer=not args.input, filename=args.input)
    pamseqobj = pam.pamSeqExp(multikmerdict=data)
    kmdkeys = pamseqobj.multikmerdict.keys()

    if args.cutoff:
        cutoff = {int(key): value for key, value in args.cutoff.items()}
        assert set(cutoff.keys()) == set(
            kmdkeys
        ), "A length listed in the --cutoff input does not match the lengths available."
    else:
        cutoff = {
            lenval: pamseqobj.find_breakpoint(length=lenval) for lenval in kmdkeys
        }
        logging.info(f"Computed cutoff: {cutoff}")
    pd_out = _create_directories(args.predict_out, map(str, kmdkeys))

    for lenval in kmdkeys:
        _mk_outputs(
            pamseqobj=pamseqobj, outdir=pd_out, lenval=lenval, cutoff=cutoff[lenval]
        )


def main(args: argparse.Namespace = None) -> None:
    """
    Executes the main workflow for identifying PAM/TAM sequences.

    This function initializes the argument parser, sets up logging, and
    executes the appropriate subcommand ('process' or 'predict') based on
    the provided command-line arguments.

    Args:
        args (argparse.Namespace, optional): Command-line arguments for
        unit testing. If not provided, arguments are parsed from sys.argv.

    Returns:
        None
    """
    parser = myparser()
    args = args or parser.parse_args()

    _logger_setup(args.log)
    logging.info("Begin processing PAM/TAM sequencing libraries")
    logging.info(args)
    if args.subcommand == "process":
        pamexpobj, run_summ = process(args)
        log_run_summary(run_summ)
        export_results(pamexpobj, args.outfile)
    elif args.subcommand == "predict":
        predict(args)


def log_run_summary(run_summ):
    """A helper function for logging FASQ processing info
    """
    logging.info(
        "Control data: {} reads processed, {} contained targets.".format(
            run_summ["cont"]["tot"], run_summ["cont"]["targets"]
        )
    )
    logging.info(
        "Experimental data: {} reads processed, {} contained targets.".format(
            run_summ["exp"]["tot"], run_summ["exp"]["targets"]
        )
    )


def export_results(pamexpobj, outfile):
    """
    Export the results from a PamiPami experiment object to an HDF5 file or standard output.

    Parameters:
        pamexpobj: An object containing the multikmerdict attribute with results to export.
        outfile (str or None): Path to the output HDF5 file. If None, writes to standard output.
    """
    logging.info("Writing results to {}".format(outfile))
    if outfile:
        tpio.export_hdf(
            multikmerdict=pamexpobj.multikmerdict,
            to_buffer=False,
            filename=outfile,
        )
    else:
        sys.stdout.buffer.write(
            tpio.export_hdf(
                multikmerdict=pamexpobj.multikmerdict,
                to_buffer=True,
                filename="tamipami.h5",
            )
        )


if __name__ == "__main__":
    main()
