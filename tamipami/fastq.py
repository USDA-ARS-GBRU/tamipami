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

from .config import config


def iterate_kmer(k: int) -> dict[str, int]:
    """
    Create a dictionary of all possible kmers of length k.
    Args:
        k: The length of the kmer.
    Returns:
        A dictionary with lexographically sorted kmers as keys and 0 as values.
    """
    bases = ["A", "C", "T", "G"]
    kmers = ["".join(p) for p in itertools.product(bases, repeat=k)]
    return {kmer: 0 for kmer in sorted(kmers)}


def merge_reads(fastq, fastq2, outfile) -> str:
    """Merge Reads and return error-corrected merged fastq

    Args:
        fastq: The path to the forward fastq file.
        fastq2: The path to reverse fastq file.
        outfile: The path to the merged fastq output.

    Returns:
        The stdout from BBmerge
    """

    try:
        parameters = [
            "bbmerge.sh",
            "in=" + fastq,
            "in2=" + fastq2,
            "out=" + outfile,
        ]
        parameters.extend(config["bbmerge"])
        print(parameters)
        p3 = subprocess.run(parameters, stderr=subprocess.PIPE)
        print(p3.stderr.decode("utf-8"))
        return p3.stderr.decode("utf-8")
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
    guide_detections = 0
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
                if pamseq in kmer_dict:
                    kmer_dict[pamseq] += 1

    return kmer_dict, tot_reads, guide_detections


def process(
    fastq: str, fastq2: str, pamlen: int, mergedfile: str, spacer: str, orientation: str
) -> dict[str, int]:
    """A function to merge the reads and count the TAM/PAM sequences

    Args:
        fastq (str):
        fastq2 (str):
        pamlen (int:
        tempdir:
        spacer:

    Returns:
    A  Dictionary containing counts of every  PAM/TAM in a sample
    """
    # mergedfile = os.path.join(tempdir, "merged.fastq.gz")
    logging.info("Merging reads and correcting sequencing errors.")
    stdout = merge_reads(fastq=fastq, fastq2=fastq2, outfile=mergedfile)
    # logging.info(stdout)
    print(stdout)
    logging.info("Counting PAM/TAM sites.")
    pamcount, tot_reads, target_detections = count_pam(
        pamlen=pamlen, spacer=spacer, fastq=mergedfile, orientation=orientation
    )
    # ref_n = check_N(list(pamcount.values()))
    # logging.info("Poisson noise is expected to be {:5.1f} % of total at N={}".format(ref_n * 100, pamlen))
    return pamcount, tot_reads, target_detections
