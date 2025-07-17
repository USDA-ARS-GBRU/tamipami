import pytest
from tamipami.fastq import iterate_kmer, merge_reads_stream, count_pam_stream, process
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io

def test_iterate_kmer_generates_all_kmers():
    k = 2
    result = iterate_kmer(k)
    expected_kmers = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    assert set(result.keys()) == set(expected_kmers)
    assert all(v == 0 for v in result.values())


def test_count_pam_stream_counts_sequences(mocker):
    # Prepare a fake FASTQ stream with a record containing the spacer
    spacer = "GATTACA"
    pamlen = 3
    orientation = "5prime"
    # The PAM is 'AAA', then the spacer
    seq = "AAAGATTACA"
    record = SeqRecord(Seq(seq), id="test", description="")
    mock_parse = mocker.patch("tamipami.fastq.SeqIO.parse", return_value=[record])
    fastq_stream = io.StringIO("irrelevant")
    result, tot_reads, guide_detections = count_pam_stream(spacer, pamlen, orientation, fastq_stream)
    assert result["AAA"] == 1
    assert guide_detections == 1
    assert tot_reads == 0  # Only one record, so enumerate starts at 0

def test_iterate_kmer_invalid_k_raises():
    with pytest.raises(ValueError):
        iterate_kmer(0)
    with pytest.raises(ValueError):
        iterate_kmer(-1)
    with pytest.raises(ValueError):
        iterate_kmer("foo")

def test_count_pam_stream_invalid_orientation_raises():
    fastq_stream = io.StringIO("irrelevant")
    with pytest.raises(ValueError):
        count_pam_stream("GATTACA", 3, "invalid_orientation", fastq_stream)

def test_process_and_count_pam_stream_error_handling(mocker):
    # Simulate merge_reads_stream raising an exception
    mocker.patch("tamipami.fastq.merge_reads_stream", side_effect=RuntimeError("merge failed"))
    with pytest.raises(RuntimeError):
        process("fq1", "fq2", 3, "GATTACA", "5prime")
    # Simulate count_pam_stream raising an exception
    mock_merge = mocker.patch("tamipami.fastq.merge_reads_stream")
    mock_proc = mocker.Mock()
    mock_proc.stdout = io.StringIO("irrelevant")
    mock_merge.return_value = mock_proc
    mocker.patch("tamipami.fastq.count_pam_stream", side_effect=RuntimeError("count failed"))
    with pytest.raises(RuntimeError):
        process("fq1", "fq2", 3, "GATTACA", "5prime")