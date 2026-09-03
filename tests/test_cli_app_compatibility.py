import os
import tempfile
import argparse
import pytest
import pandas as pd
import numpy as np

from tamipami import pam, tpio, cli, cli_utils


def create_sample_pamseqexp():
    """Helper to generate a dummy pamSeqExp object with 3 lengths (3, 4, 5)."""
    # Create simple mock kmer counts for ctl and exp (length 5)
    bases = ["A", "C", "G", "T"]
    import itertools
    kmers_5 = ["".join(p) for p in itertools.product(bases, repeat=5)]

    np.random.seed(42)
    ctl_counts = {k: int(np.random.randint(50, 100)) for k in kmers_5}
    exp_counts = {k: int(np.random.randint(10, 100)) for k in kmers_5}
    # Deplete some kmers
    for k in kmers_5[:20]:
        exp_counts[k] = 0

    return pam.pamSeqExp(ctl=ctl_counts, exp=exp_counts, position="3prime")


def test_cli_length_arg_default_and_choices():
    """Test CLI parser defaults and valid/invalid --length choices."""
    parser = cli.myparser()

    # Default length should be 5
    args = parser.parse_args(["process", "-c", "c1.fq", "-e", "e1.fq"])
    assert args.length == 5

    # Length 3, 4, 5 should be accepted
    for l in [3, 4, 5]:
        args = parser.parse_args(["process", "-c", "c1.fq", "-e", "e1.fq", "--length", str(l)])
        assert args.length == l

    # Length 6 should be rejected (out of 3-5 range)
    with pytest.raises(SystemExit):
        parser.parse_args(["process", "-c", "c1.fq", "-e", "e1.fq", "--length", "6"])


def test_cli_cutoff_validator():
    """Test --cutoff validator accepts keys 3, 4, 5 and rejects out-of-range keys."""
    valid_json = '{"3": 0.5, "4": 1.0, "5": 1.5}'
    res = cli_utils.cutoff_arg_validator(valid_json)
    assert res == {3: 0.5, 4: 1.0, 5: 1.5}

    invalid_json_key = '{"6": 0.5}'
    with pytest.raises(argparse.ArgumentTypeError):
        cli_utils.cutoff_arg_validator(invalid_json_key)


def test_cli_mk_outputs_and_logo_filtering(mocker):
    """Test _mk_outputs generates all expected files and filters kmers >= cutoff for logo."""
    pamexpobj = create_sample_pamseqexp()

    with tempfile.TemporaryDirectory() as tmpdir:
        out_subdir = os.path.join(tmpdir, "3")
        os.makedirs(out_subdir, exist_ok=True)

        # Mock make_logo to verify parameters
        spy_make_logo = mocker.spy(pamexpobj, "make_logo")

        cutoff_val = 1.0
        cli._mk_outputs(pamseqobj=pamexpobj, outdir=tmpdir, lenval=3, cutoff=cutoff_val)

        # Check logo was called with above=False and score_type="zscore"
        spy_make_logo.assert_called_once()
        call_kwargs = spy_make_logo.call_args.kwargs
        assert call_kwargs["cutoff"] == cutoff_val
        assert call_kwargs["above"] is False
        assert call_kwargs["score_type"] == "zscore"

        # Check output files exist
        assert os.path.exists(os.path.join(out_subdir, "data.csv"))
        assert os.path.exists(os.path.join(out_subdir, "histogram.html"))
        assert os.path.exists(os.path.join(out_subdir, "degenerate_pam_tam.csv"))

        # Verify degenerate output file content
        degen_df = pd.read_csv(os.path.join(out_subdir, "degenerate_pam_tam.csv"))
        assert "PAM/TAM site" in degen_df.columns


def test_app_cli_hdf5_bidirectional_compatibility():
    """Test HDF5 exported by App/CLI can be re-imported and predicted by CLI."""
    pamexpobj = create_sample_pamseqexp()

    # 1. Export HDF5 buffer as App does
    hdf_bytes = tpio.export_hdf(pamexpobj.multikmerdict, to_buffer=True, filename="tamipami.h5")
    assert isinstance(hdf_bytes, bytes)

    # 2. Save HDF5 buffer to disk
    with tempfile.TemporaryDirectory() as tmpdir:
        h5_path = os.path.join(tmpdir, "app_exported.h5")
        with open(h5_path, "wb") as f:
            f.write(hdf_bytes)

        # 3. Import HDF5 via CLI tpio.import_hdf
        imported_dict = tpio.import_hdf(from_buffer=False, filename=h5_path)
        assert set(imported_dict.keys()) == {3, 4, 5}

        # Verify column consistency across all kmer lengths
        expected_cols = [
            "kmers", "ctl_raw", "exp_raw", "ctl_clr", "exp_clr",
            "diff", "zscore", "pvalue", "p_adjust_BH",
            "is_statistically_cut", "true_percent_depleted", "shrunk_percent_depleted"
        ]
        for klen in [3, 4, 5]:
            for col in expected_cols:
                assert col in imported_dict[klen].columns

        # 4. Run CLI predict command using the App-exported HDF5 file
        predict_outdir = os.path.join(tmpdir, "predict_results")
        parser = cli.myparser()
        args = parser.parse_args([
            "predict",
            "--input", h5_path,
            "--cutoff", '{"3": 0.5, "4": 0.5, "5": 0.5}',
            "--predict_out", predict_outdir
        ])

        cli.predict(args)

        # Check outputs exist for each length
        for klen in ["3", "4", "5"]:
            kdir = os.path.join(predict_outdir, klen)
            assert os.path.exists(os.path.join(kdir, "data.csv"))
            assert os.path.exists(os.path.join(kdir, "histogram.html"))
            assert os.path.exists(os.path.join(kdir, "degenerate_pam_tam.csv"))
