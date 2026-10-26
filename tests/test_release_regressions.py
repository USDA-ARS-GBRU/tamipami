import itertools
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

from tamipami import candidates, cli, degenerate, fastq, pam, tpio


def _kmer_counts(lam, depleted_frac, seed):
    kmers = ["".join(p) for p in itertools.product("ACGT", repeat=5)]
    rng = np.random.default_rng(seed)
    ctl = {k: int(v) for k, v in zip(kmers, rng.poisson(lam, len(kmers)))}
    exp = {
        k: (int(v) if i % 7 else int(v * depleted_frac))
        for i, (k, v) in enumerate(zip(kmers, rng.poisson(lam, len(kmers))))
    }
    return ctl, exp


def test_low_depth_strong_depletion_gives_finite_zscores():
    ctl, exp = _kmer_counts(lam=3, depleted_frac=0.0, seed=1)
    obj = pam.pamSeqExp(ctl=ctl, exp=exp, position="3prime")
    for df in obj.multikmerdict.values():
        assert np.isfinite(df["zscore"]).all()


def test_identical_inputs_give_identical_results():
    ctl, exp = _kmer_counts(lam=200, depleted_frac=0.01, seed=0)
    a = pam.pamSeqExp(ctl=ctl, exp=exp, position="5prime").multikmerdict[5]
    b = pam.pamSeqExp(ctl=ctl, exp=exp, position="5prime").multikmerdict[5]
    pd.testing.assert_frame_equal(a, b)


def test_export_hdf_overwrites_existing_file():
    df = pd.DataFrame({"kmers": ["AAA", "CCC"], "zscore": [1.0, 2.0]})
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "out.h5")
        tpio.export_hdf({3: df}, to_buffer=False, filename=path)
        tpio.export_hdf({3: df}, to_buffer=False, filename=path)
        assert len(tpio.import_hdf(filename=path)[3]) == 2


def test_export_hdf_buffer_ignores_file_on_disk(monkeypatch):
    df = pd.DataFrame({"kmers": ["AAA", "CCC"], "zscore": [1.0, 2.0]})
    with tempfile.TemporaryDirectory() as tmp:
        monkeypatch.chdir(tmp)
        tpio.export_hdf({3: df}, to_buffer=False, filename="tamipami.h5")
        data = tpio.export_hdf({3: df}, to_buffer=True, filename="tamipami.h5")
        with open("buf.h5", "wb") as f:
            f.write(data)
        assert len(tpio.import_hdf(filename="buf.h5")[3]) == 2
        # The on-disk file is left untouched
        assert len(tpio.import_hdf(filename="tamipami.h5")[3]) == 2


def test_dense_kmer_set_solves_quickly():
    seqs = ["".join(p) for p in itertools.product("ACGT", repeat=5)]
    assert degenerate.seqs_to_degenerates(seqs) == ["NNNNN"]
    with pytest.raises(candidates.CandidateLimitExceeded):
        candidates.positional_subset_candidates(seqs, max_combinations=1000)


def test_positional_candidates_keep_largest_when_capped():
    seqs = ["".join(p) for p in itertools.product("ACGT", repeat=3)]
    kept = candidates.positional_subset_candidates(seqs, max_candidates=10)
    assert len(kept) == 10
    assert ("NNN", set(seqs)) in kept


def test_solve_never_worse_than_tree_cover():
    rng = np.random.default_rng(3)
    allk = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
    seqs = sorted(rng.choice(allk, size=150, replace=False).tolist())
    tree = candidates.tree_induced_candidates(seqs)
    tree_codes = degenerate.minimal_exact_cover(seqs, tree, time_limit=5)
    codes = degenerate.seqs_to_degenerates(seqs, deterministic_time=2)
    assert degenerate._cover_score(codes) <= degenerate._cover_score(tree_codes)
    covered = sorted(s for c in codes for s in candidates.expand_iupac(c))
    assert covered == seqs  # exact cover: every k-mer exactly once


def test_serve_uses_packaged_app_path(mocker):
    run = mocker.patch("subprocess.run")
    mocker.patch("tamipami.cli.logger_setup")
    args = cli.myparser().parse_args(["serve"])
    cli.main(args)
    script = run.call_args.args[0][4]
    assert os.path.isabs(script)
    assert os.path.exists(script)


def test_bbmerge_failure_raises(mocker, tmp_path):
    fq = tmp_path / "reads.fq"
    fq.write_text("")
    mocker.patch(
        "tamipami.fastq.merge_reads_stream",
        side_effect=lambda fastq, fastq2, stderr: subprocess.Popen(
            [sys.executable, "-c", "import sys; sys.stderr.write('boom'); sys.exit(3)"],
            stdout=subprocess.PIPE,
            stderr=stderr,
            text=True,
        ),
    )
    with pytest.raises(RuntimeError, match="exit code 3"):
        fastq.process(str(fq), str(fq), 3, "GATTACA", "5prime")


@pytest.mark.skipif(shutil.which("bbmerge.sh") is None, reason="bbmerge.sh not installed")
def test_paired_end_paths_with_shell_metacharacters(tmp_path):
    src = os.path.join(os.path.dirname(__file__), os.pardir, "data", "fastq")
    d = tmp_path / "Example (2); echo pwned"
    d.mkdir()
    for r in ("R1", "R2"):
        shutil.copy(os.path.join(src, f"sample1k_LibraryRTW574_{r}_003.fastq.gz"), d / f"{r}.fastq.gz")
    _, tot, targets = fastq.process(
        str(d / "R1.fastq.gz"), str(d / "R2.fastq.gz"), 3, "CTGATGGTCCATGTCTGTTACTC", "5prime"
    )
    assert tot > 0 and targets > 0


def test_paired_end_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        fastq.process(str(tmp_path / "a.fq"), str(tmp_path / "b.fq"), 3, "GATTACA", "5prime")


def test_degenerate_solve_is_reproducible():
    rng = np.random.default_rng(4)
    allk = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
    seqs = sorted(rng.choice(allk, size=150, replace=False).tolist())
    first = degenerate.seqs_to_degenerates(seqs, deterministic_time=5)
    assert all(degenerate.seqs_to_degenerates(seqs, deterministic_time=5) == first for _ in range(2))
