import os

# -*- coding: utf-8 -*-
"""
fastq: A tamipami module for processing FASTQs from  PAM/TAM sequencing libraries

"""
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

import subprocess
import logging
import re
from Bio import SeqIO

from tamipami.config import config

def merge_reads_stream(fastq: str, fastq2: str) -> subprocess.Popen:
    """
    Merge paired-end FASTQ files using BBmerge and stream merged reads to stdout.

    Args:
        fastq: The path to the forward FASTQ file.
        fastq2: The path to the reverse FASTQ file.

    Returns:
        A subprocess.Popen object with stdout as a stream of merged FASTQ records.
    """
    parameters = [
        "bbmerge.sh",
        "in=" + fastq,
        "in2=" + fastq2,
        "out=stdout.fastq",# Stream merged reads to stdout
    ]
    parameters.extend(config["bbmerge"])
    logging.info(parameters)
    proc = subprocess.Popen(parameters, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return proc

def count_pam_stream(
    spacer: str, pamlen: int, orientation: str, fastq_stream
) -> tuple[dict[str, int], int, int]:
    """
    Count PAM/TAM sequences in a streamed FASTQ file-like object.

    Args:
        spacer: The DNA sequence matching the guide sequence.
        pamlen: The length of the PAM/TAM sequences.
        orientation: '5prime' (5'-PAM-Spacer-3') or '3prime' (5'-Spacer-PAM-3').
        fastq_stream: File-like object streaming FASTQ records.

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
    tot_reads = 0
    try:
        for tot_reads, record in enumerate(SeqIO.parse(fastq_stream, "fastq")):
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
    except Exception as e:
        logging.error("Error during streaming PAM/TAM counting: %s", e)
        raise

def process(
    fastq: str, fastq2: str, pamlen: int, spacer: str, orientation: str
) -> tuple[dict[str, int], int, int]:
    """
    Merge reads from FASTQ files and count PAM/TAM sequences on the fly.

    Args:
        fastq (str): Path to the forward FASTQ file.
        fastq2 (str): Path to the reverse FASTQ file.
        pamlen (int): Length of the PAM/TAM sequences.
        spacer (str): DNA sequence matching the guide sequence.
        orientation (str): Orientation of the PAM/TAM sequence, either '5prime' or '3prime'.

    Returns:
        dict[str, int]: Dictionary containing counts of each PAM/TAM sequence in the sample.
    """
    if orientation not in {"5prime", "3prime"}:
        raise ValueError("`orientation` must be '5prime' or '3prime'.")
    logging.info("Merging reads and streaming to PAM/TAM counter.")
    try:
        proc = merge_reads_stream(fastq=fastq, fastq2=fastq2)
        pamcount, tot_reads, target_detections = count_pam_stream(
            pamlen=pamlen, spacer=spacer, orientation=orientation, fastq_stream=proc.stdout
        )
        proc.stdout.close()
        proc.wait()
    except Exception as e:
        logging.error(f"Error during streaming merge/count: {e}")
        raise
    return pamcount, tot_reads, target_detections