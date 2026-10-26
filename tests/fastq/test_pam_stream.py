import pytest
import io

from tamipami.fastq import count_pam_stream

def make_fastq(records):
    """Helper to create a FASTQ string from a list of (name, seq, qual) tuples."""
    return "\n".join(
        f"@{name}\n{seq}\n+\n{qual}" for name, seq, qual in records
    ) + "\n"

def test_count_pam_stream_5prime_correct_count():
    # Guide: GGG, PAM length: 2, orientation: 5prime
    # Sequence: XXPAMGGGYY, PAM is at positions 2:4, guide at 4:7
    #records = [
    #    ("r1", "AACCGGGTT", "IIIIIIIII"),  # PAM: CC, guide: GGG
    #    ("r2", "TTCCGGGAA", "IIIIIIIII"),  # PAM: CC, guide: GGG
    #    ("r3", "AACCGGGTT", "IIIIIIIII"),  # PAM: CC, guide: GGG
    #    ("r4", "AAGGGTTCC", "IIIIIIIII"),  # PAM: AA, guide: GGG
    #]
    records = [
        ("r1", 'AACCCGGTT', "IIIIIIIII"),  # PAM: CC, guide: GGG
        ("r2", 'TTCCCGGAA', "IIIIIIIII"),  # PAM: CC, guide: GGG
        ("r3", 'AACCCGGTT', "IIIIIIIII"),  # PAM: CC, guide: GGG
        ("r4", 'GGAACCCTT', "IIIIIIIII"),  # PAM: AA, guide: GGG
    ]
    fastq = io.StringIO(make_fastq(records))
    kmer_counts, tot_reads, guide_detections = count_pam_stream("GGG", 2, "5prime", fastq)
    assert tot_reads == 4
    assert guide_detections == 4
    # Only CC and AA should be counted
    assert kmer_counts["CC"] == 3
    assert kmer_counts["AA"] == 1
    # All other k-mers should be zero
    for kmer, count in kmer_counts.items():
        if kmer not in {"CC", "AA"}:
            assert count == 0

def test_count_pam_stream_3prime_correct_count():
    # Guide: GGG, PAM length: 2, orientation: 3prime
    # Sequence: XXGGGPAMYY, PAM is at positions 5:7, guide at 2:5
    records = [
        ("r1", 'TTGGCCCTT', "IIIIIIIII"),  # PAM: CC, guide: GGG
        ("r2", 'TTAACCCAA', "IIIIIIIII"),  # PAM: TT, guide: GGG
        ("r3", 'TTGGCCCTT', "IIIIIIIII"),  # PAM: CC, guide: GGG
        ("r4", 'TTCCCCCTT', "IIIIIIIII"),  # PAM: GG, guide: GGG
    ]
    fastq = io.StringIO(make_fastq(records))
    kmer_counts, tot_reads, guide_detections = count_pam_stream("GGG", 2, "3prime", fastq)
    assert tot_reads == 4
    assert guide_detections == 4
    assert kmer_counts["CC"] == 2
    assert kmer_counts["TT"] == 1
    assert kmer_counts["GG"] == 1
    for kmer, count in kmer_counts.items():
        if kmer not in {"CC", "TT", "GG"}:
            assert count == 0

def test_count_pam_stream_no_guide_sequence():
    # No guide sequence present
    records = [
        ("r1", "AACCTTGG", "IIIIIIII"),
        ("r2", "TTCCAATT", "IIIIIIII"),
    ]
    fastq = io.StringIO(make_fastq(records))
    kmer_counts, tot_reads, guide_detections = count_pam_stream("GGG", 2, "5prime", fastq)
    assert tot_reads == 2
    assert guide_detections == 0
    assert all(count == 0 for count in kmer_counts.values())

def test_count_pam_stream_invalid_orientation():
    records = [
        ("r1", "AACCGGGTT", "IIIIIIIII"),
    ]
    fastq = io.StringIO(make_fastq(records))
    with pytest.raises(ValueError, match="`orientation` must be '5prime' or '3prime'"):
        count_pam_stream("GGG", 2, "invalid", fastq)

def test_count_pam_stream_malformed_fastq():
    # Malformed FASTQ: missing quality line
    malformed_fastq = "@r1\nAACCGGGTT\n+\nIIIIIIIII\n@r2\nTTCCAATT\n+\n"  # missing quality for r2
    fastq = io.StringIO(malformed_fastq)
    with pytest.raises(Exception):
        count_pam_stream("GGG", 2, "5prime", fastq)


@pytest.mark.parametrize("spacer", ["", "(A+)+C", "GGN", "GG G", None])
def test_count_pam_stream_rejects_invalid_spacer(spacer):
    fastq = io.StringIO(make_fastq([("r1", "AACCCGGTT", "IIIIIIIII")]))
    with pytest.raises(ValueError):
        count_pam_stream(spacer, 2, "5prime", fastq)


def test_count_pam_stream_normalizes_spacer_case_and_whitespace():
    fastq = io.StringIO(make_fastq([("r1", "AACCCGGTT", "IIIIIIIII")]))
    kmer_counts, _, guide_detections = count_pam_stream(" ggg ", 2, "5prime", fastq)
    assert guide_detections == 1
    assert kmer_counts["CC"] == 1


def test_count_pam_stream_5prime_truncated_pam_not_counted():
    # Reverse complement is GGGTTCCAA: the spacer is at the read start, so no PAM precedes it
    fastq = io.StringIO(make_fastq([("r1", "TTGGAACCC", "IIIIIIIII")]))
    kmer_counts, tot_reads, guide_detections = count_pam_stream("GGG", 2, "5prime", fastq)
    assert tot_reads == 1
    assert guide_detections == 1
    assert sum(kmer_counts.values()) == 0
