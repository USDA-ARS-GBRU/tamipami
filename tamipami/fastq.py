# -*- coding: utf-8 -*-
"""
fastq: A tamipami module for processing FASTQs from  PAM/TAM sequencing libraries

"""
import os
import subprocess
import tempfile
import logging
import re
import itertools
import gzip
from typing import Optional

from Bio import SeqIO

from tamipami.config import config

logger = logging.getLogger(__name__)

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


def validate_spacer(spacer: str) -> str:
    """
    Normalize a spacer sequence and verify it contains only A, C, G and T.

    Args:
        spacer: The spacer sequence supplied by the user.

    Returns:
        The upper-cased spacer with surrounding whitespace removed.

    Raises:
        ValueError: If the spacer is empty or contains characters other than A, C, G, T.
    """
    if not isinstance(spacer, str):
        raise ValueError("A spacer sequence must be provided, or a library selected.")
    spacer = spacer.strip().upper()
    if not re.fullmatch(r"[ACGT]+", spacer):
        raise ValueError(
            "The spacer sequence must be non-empty and contain only A, C, G and T."
        )
    return spacer


def _safe_link(path: str, linkdir: str, name: str) -> str:
    """
    Symlink `path` into `linkdir` under a shell-safe name, keeping a .gz suffix.

    bbmerge.sh passes its arguments through `eval`, so paths containing spaces,
    parentheses or other shell metacharacters must never reach it directly.
    """
    linkname = os.path.join(linkdir, name + (".fastq.gz" if path.endswith(".gz") else ".fastq"))
    os.symlink(os.path.abspath(path), linkname)
    return linkname


def merge_reads_stream(fastq: str, fastq2: str, stderr=subprocess.PIPE) -> subprocess.Popen:
    """
    Merge paired-end FASTQ files using BBmerge and stream merged reads to stdout.

    Args:
        fastq: The path to the forward FASTQ file.
        fastq2: The path to the reverse FASTQ file.
        stderr: Where BBmerge's stderr is sent. Pass a file object to avoid filling the
            pipe buffer while stdout is being consumed.

    Returns:
        A subprocess.Popen object with stdout as a stream of merged FASTQ records.
    """
    parameters = [
        "bbmerge.sh",
        "in1=" + fastq,
        "in2=" + fastq2,
        "out=stdout.fastq",  # Stream merged reads to stdout
        "interleaved=f",
    ]
    parameters.extend(config["bbmerge"])
    logger.info(parameters)
    proc = subprocess.Popen(
        parameters, stdout=subprocess.PIPE, stderr=stderr, text=True
    )
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
    spacer = validate_spacer(spacer)
    pamlen = int(pamlen)

    kmer_dict = iterate_kmer(pamlen)
    guide_detections = 0
    tot_reads = 0
    try:
        seqrecords = SeqIO.parse(fastq_stream, "fastq")
        for record in seqrecords:
            tot_reads += 1
            seqstr = str(record.reverse_complement().seq)
            spacerstart = seqstr.find(spacer)
            if spacerstart == -1:
                continue
            guide_detections += 1
            if orientation == "5prime":
                pamstart = spacerstart - pamlen
                if pamstart < 0:
                    continue  # PAM truncated by the start of the read
                pamseq = seqstr[pamstart:spacerstart]
            else:
                spacerend = spacerstart + len(spacer)
                pamseq = seqstr[spacerend : spacerend + pamlen]
            if pamseq in kmer_dict:
                kmer_dict[pamseq] += 1
        return kmer_dict, tot_reads, guide_detections
    except Exception as e:
        logger.error("Error during streaming PAM/TAM counting: %s", e)
        raise


def process(
    fastq: str, fastq2: Optional[str], pamlen: int, spacer: str, orientation: str
) -> tuple[dict[str, int], int, int]:
    """
    Merge reads from paired-end FASTQ files (via BBmerge) or stream single-end FASTQ directly,
    then count PAM/TAM sequences on the fly.

    Args:
        fastq (str): Path to the forward FASTQ file.
        fastq2 (Optional[str]): Path to the reverse FASTQ file for paired-end runs. If None, single-end is assumed.
        pamlen (int): Length of the PAM/TAM sequences.
        spacer (str): DNA sequence matching the guide sequence.
        orientation (str): Orientation of the PAM/TAM sequence, either '5prime' or '3prime'.

    Returns:
        tuple[dict[str, int], int, int]: A tuple of (kmer counts, total reads processed, guide detections).
    """
    if orientation not in {"5prime", "3prime"}:
        raise ValueError("`orientation` must be '5prime' or '3prime'.")
    spacer = validate_spacer(spacer)

    try:
        if fastq2:
            logger.info(
                "Paired-end input detected; merging reads and streaming to PAM/TAM counter."
            )
            for path in (fastq, fastq2):
                if not os.path.isfile(path):
                    raise FileNotFoundError(f"Input file not found: {path}")
            # BBmerge stderr goes to a temp file so a chatty run cannot fill the pipe and deadlock
            with tempfile.TemporaryDirectory(prefix="tamipami_bbmerge_") as linkdir, \
                    tempfile.TemporaryFile(mode="w+") as errfile:
                in1 = _safe_link(fastq, linkdir, "in1")
                in2 = _safe_link(fastq2, linkdir, "in2")
                proc = merge_reads_stream(fastq=in1, fastq2=in2, stderr=errfile)
                try:
                    pamcount, tot_reads, target_detections = count_pam_stream(
                        pamlen=pamlen,
                        spacer=spacer,
                        orientation=orientation,
                        fastq_stream=proc.stdout,
                    )
                finally:
                    proc.stdout.close()
                    returncode = proc.wait()
                    errfile.seek(0)
                    stderr_output = errfile.read()
                if returncode != 0:
                    logger.error(f"BBmerge error: {stderr_output}")
                    raise RuntimeError(
                        f"BBmerge failed with exit code {returncode}. Please verify your input files."
                    )
                if stderr_output:
                    logger.info(f"BBmerge output: {stderr_output}")

        else:
            logger.info(
                "Single-end input detected; streaming reads directly to PAM/TAM counter."
            )
            if fastq.endswith(".gz"):
                fh = gzip.open(fastq, "rt")
            else:
                fh = open(fastq, "rt")
            try:
                pamcount, tot_reads, target_detections = count_pam_stream(
                    pamlen=pamlen,
                    spacer=spacer,
                    orientation=orientation,
                    fastq_stream=fh,
                )
            finally:
                fh.close()
    except Exception as e:
        logger.error(f"Error during streaming merge/count: {e}")
        raise
    return pamcount, tot_reads, target_detections
