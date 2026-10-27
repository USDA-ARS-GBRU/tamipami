#!/usr/bin/env python
"""Degenerate-solver time, memory and result quality by selection size.

For k = 3, 4, 5, the k-mers in a `tamipami process` HDF5 file are ranked by z-score
and the top fraction (or the automatic cutoff) is passed to
tamipami.degenerate.seqs_to_degenerates with its default settings, as the app and
CLI call it. Each solve runs in a separate process so peak memory is measured per
solve and runaway solves can be killed.

--pythonpath runs the solve against another checkout of tamipami, e.g. a git
worktree of an older commit, to compare versions.

Example:
  python run_selection_sweep.py --hdf results/example_kmers.h5 --out results/selection_current.jsonl
  git worktree add /tmp/tp_old 101284e
  python run_selection_sweep.py --hdf results/example_kmers.h5 --pythonpath /tmp/tp_old \\
      --out results/selection_101284e.jsonl
"""
import argparse
import json
import os
import subprocess
import sys

FRACTIONS = ("auto", 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)


def worker(hdf, k, fraction):
    """Solve one selection and print a JSON line (runs in a child process)."""
    import logging
    import resource
    import time
    import warnings

    warnings.filterwarnings("ignore")
    from tamipami import degenerate, pam, tpio

    statuses, limited = [], []

    class H(logging.Handler):
        def emit(self, record):
            msg = record.getMessage()
            if "status=" in msg:
                statuses.append(msg.split("status=")[1].split()[0].replace("CpSolverStatus.", ""))
            if "Skipping positional" in msg or "Kept the" in msg:
                limited.append(msg)

    for name in ("tamipami.degenerate", "tamipami.candidates"):
        logging.getLogger(name).addHandler(H())
        logging.getLogger(name).setLevel(logging.INFO)

    data = tpio.import_hdf(filename=hdf)
    df = data[k].sort_values("zscore", ascending=False)
    if fraction == "auto":
        cutoff = pam.pamSeqExp(multikmerdict=data).find_breakpoint(k)
        seqs = sorted(df[df.zscore >= cutoff].kmers)
    else:
        seqs = sorted(df.kmers.iloc[: max(1, round(float(fraction) * len(df)))])

    t0 = time.time()
    codes = degenerate.seqs_to_degenerates(seqs)
    elapsed = time.time() - t0
    iupac_size = {c: len(v) for c, v in degenerate.IUPAC_CODES.items()}
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_mb = rss / 2**20 if sys.platform == "darwin" else rss / 2**10
    print(json.dumps(dict(
        k=k, fraction=fraction, n=len(seqs), seconds=round(elapsed, 2), peak_rss_mb=round(rss_mb),
        stage_statuses=statuses, candidates_limited=bool(limited),
        patterns=len(codes),
        degenerate_positions=sum(iupac_size[ch] > 1 for c in codes for ch in c),
        bases=sum(iupac_size[ch] for c in codes for ch in c),
        codes=codes if len(codes) <= 10 else codes[:10] + ["..."],
    )))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hdf", required=True, help="tamipami process HDF5 output")
    p.add_argument("--out", required=True, help="JSON-lines output file (overwritten)")
    p.add_argument("--pythonpath", help="tamipami checkout to import instead of the installed package")
    p.add_argument("--timeout", type=float, default=300.0, help="seconds before a solve is killed")
    p.add_argument("--label", default="", help="version label stored with each result")
    p.add_argument("--worker", nargs=3, metavar=("HDF", "K", "FRACTION"), help=argparse.SUPPRESS)
    args = p.parse_args()

    env = dict(os.environ)
    if args.pythonpath:
        env["PYTHONPATH"] = os.path.abspath(args.pythonpath)
    hdf = os.path.abspath(args.hdf)
    with open(args.out, "w") as out:
        for k in (3, 4, 5):
            for fraction in FRACTIONS:
                cmd = [sys.executable, os.path.abspath(__file__), "--hdf", hdf, "--out", os.devnull,
                       "--worker", hdf, str(k), str(fraction)]
                try:
                    res = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=args.timeout)
                    lines = [ln for ln in res.stdout.splitlines() if ln.startswith("{")]
                    rec = json.loads(lines[-1]) if lines else dict(
                        k=k, fraction=fraction, error=(res.stderr.strip().splitlines() or ["no output"])[-1])
                except subprocess.TimeoutExpired:
                    rec = dict(k=k, fraction=fraction, error=f"killed after {args.timeout:.0f} s")
                rec["version"] = args.label
                out.write(json.dumps(rec) + "\n")
                out.flush()
                print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--worker" in args:
        i = args.index("--worker")
        hdf, k, fraction = args[i + 1 : i + 4]
        worker(hdf, int(k), fraction if fraction == "auto" else float(fraction))
    else:
        main()
