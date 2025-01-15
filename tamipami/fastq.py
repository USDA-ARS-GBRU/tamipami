import os

# -*- coding: utf-8 -*-
"""
fastq: A tamipami module for processing FASTQs from  PAM/TAM sequencing libraries

"""
import itertools
import subprocess
import logging
import re
import itertools
import gzip

from Bio import SeqIO

from tamipami.config import config


def iterate_kmer(k: int) -> dict[str, int]:
    """
    Generate a dictionary of all possible k-mers of a specified length.

    Args:
        k (int): The length of the k-mer, must be a positive integer.

    Returns:
        dict[str, int]: A dictionary with lexicographically sorted k-mers as keys and 0 as values.

    Raises:
        ValueError: If `k` is not a positive integer.
    """
    if not isinstance(k, int) or k < 1:
        raise ValueError("`k` must be a positive integer.")
    bases = ["A", "C", "T", "G"]
    kmers = ["".join(p) for p in itertools.product(bases, repeat=k)]
    return {kmer: 0 for kmer in sorted(kmers)}


def merge_reads(fastq: str, fastq2: str, outfile: str) -> str:
    """
    Merge paired-end FASTQ files using BBmerge and return the error-corrected, merged FASTQ file.

    This function utilizes the BBmerge tool to merge forward and reverse FASTQ files
    into a single output file. It logs the command parameters and captures both
    standard output and error output from the BBmerge process.

    Args:
        fastq: The path to the forward FASTQ file.
        fastq2: The path to the reverse FASTQ file.
        outfile: The path to the output merged FASTQ file.

    Returns:
        The standard output from the BBmerge process as a string. Returns an empty
        string if an error occurs during the merging process.
    """

    fastq = os.path.join(fastq)
    fastq2 = os.path.join(fastq2)
    outfile = os.path.join(outfile)

    try:
        parameters = [
            "bbmerge.sh",
            "in=" + fastq,
            "in2=" + fastq2,
            "out=" + outfile,
        ]
        parameters.extend(config["bbmerge"])
        logging.info(parameters)
        p3 = subprocess.run(
            parameters, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        stdout_output = p3.stdout.decode("utf-8")
        stderr_output = p3.stderr.decode("utf-8")
        logging.info(stderr_output)
        return stdout_output
    except subprocess.CalledProcessError as e:
        logging.error(f"Subprocess error: {e.stderr.decode('utf-8')}")
        raise e
    except RuntimeError as e:
        logging.warning("could not perform read merging with BBmerge")
        raise e



def count_pam(
    spacer: str, fastq: str, pamlen: int, orientation: str
) -> tuple[dict[str, int], int, int]:
    """
    Count PAM/TAM sequences in a FASTQ file.

    Args:
        spacer: The DNA sequence matching the guide sequence.
        fastq: The path to a merged FASTQ file.
        pamlen: The length of the PAM/TAM sequences.
        orientation: '5prime' (5'-PAM-Spacer-3') or '3prime' (5'-Spacer-PAM-3').

    Returns:
        A tuple containing:
        - A dictionary with kmer counts.
        - The total number of reads processed.
        - The number of guide sequence detections.
    """
    if orientation not in {"5prime", "3prime"}:
        raise ValueError("`orientation` must be '5prime' or '3prime'.")

    kmer_dict = iterate_kmer(pamlen)
    guide_detections = 0
    try:
        tot_reads = 0
        with gzip.open(fastq, "rt") as handle:
            for tot_reads, record in enumerate(SeqIO.parse(handle, "fastq")):
                seqstr = str(record.seq)
                result = re.search(spacer, seqstr)
                if result:
                    if orientation == "5prime":
                        spacerstart = result.start()
                        pamstart = spacerstart - int(pamlen)
                        pamseq = seqstr[pamstart:spacerstart]
                        guide_detections += 1
                    elif orientation == "3prime":
                        spacerend = result.end()
                        pamend = spacerend + int(pamlen)
                        pamseq = seqstr[spacerend:pamend]
                        guide_detections += 1
                    if pamseq and pamseq in kmer_dict:
                        kmer_dict[pamseq] += 1

        return kmer_dict, tot_reads, guide_detections
    except (ValueError, TypeError) as e:
        logging.error("Could not read the merged FASTQ file. please verify the input FASTQ and merged FASTQ files ere not empty.")
        raise e
    except (gzip.BadGzipFile, OSError, EOFError) as e:
        logging.error("Could not open the gzipped, merged FASTQ file")
        raise e
    except Exception as e:
        logging.error("could not open and count PAM sites in the Merged FASTQ file")
        raise e
    else:
        logging.error("An unexpected error occurred")
        raise e

def process(
    fastq: str, fastq2: str, pamlen: int, mergedfile: str, spacer: str, orientation: str
) -> tuple[dict[str, int], int, int]:
    """Merge reads from FASTQ files and count PAM/TAM sequences.

    Args:
        fastq (str): Path to the forward FASTQ file.
        fastq2 (str): Path to the reverse FASTQ file.
        pamlen (int): Length of the PAM/TAM sequences.
        mergedfile (str): Path to the output merged FASTQ file.
        spacer (str): DNA sequence matching the guide sequence.
        orientation (str): Orientation of the PAM/TAM sequence, either '5prime' or '3prime'.

    Returns:
        dict[str, int]: Dictionary containing counts of each PAM/TAM sequence in the sample.
    """
    if orientation not in {"5prime", "3prime"}:
        raise ValueError("`orientation` must be '5prime' or '3prime'.")
    logging.info("Merging reads and correcting sequencing errors.")
    try:
        stdout = merge_reads(fastq=fastq, fastq2=fastq2, outfile=mergedfile)
        print(stdout)
    except Exception as e:
        logging.error(f"Error during merging reads: {e}")
        return {}

    logging.info("Counting PAM/TAM sites.")
    try:
        pamcount, tot_reads, target_detections = count_pam(
            pamlen=pamlen, spacer=spacer, fastq=mergedfile, orientation=orientation
        )
    except Exception as e:
        logging.error(f"Error during PAM/TAM counting: {e}")
        return {}
    return pamcount, tot_reads, target_detections
