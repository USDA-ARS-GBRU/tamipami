#!/usr/bin/env python3
"""Benchmark the TamiPami CLI (`process` then `predict`) on each experiment, N runs each.

Versions come from git, never from this file. By default two targets are benchmarked: the most
recent release tag and the current working tree. Use ``--ref`` (repeatable) for any other tag,
branch or commit. Results go to ``results/cli_<label>.csv`` with a ``.meta.json`` sidecar.
"""
import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import time
from functools import lru_cache
from pathlib import Path

import psutil

import bench_common as bc

CUTOFF_VALUE = 2
FIELDS = ["Experiment", "Enzyme", "Run", "process_time", "predict_time", "process_max_mem",
          "predict_max_mem", "process_avg_cpu", "predict_avg_cpu", "Control_reads", "Exp_reads"]


@lru_cache(maxsize=None)
def count_fastq_pairs(path: str) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def cutoff_for(h5_path) -> str:
    """Fixed z-score cutoff for every k-mer length stored in the HDF5 file. Some releases require
    --cutoff to name every available length, so read the lengths rather than hard coding them."""
    import pandas as pd
    with pd.HDFStore(str(h5_path), mode="r") as store:
        lengths = sorted(int(k.strip("/").split("_")[1]) for k in store.keys())
    return json.dumps({str(k): CUTOFF_VALUE for k in lengths})


def run_measured(cmd, env):
    """Run cmd; return (wall_seconds, peak_rss_MB, avg_cpu_percent) over the whole process tree
    (so the bbmerge.sh/Java child is included). CPU% is total CPU time / wall time, so 100 means
    one core fully busy and values above 100 mean multiple cores."""
    start = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=env)
    root = psutil.Process(proc.pid)
    peak, cpu_by_pid = 0, {}
    while proc.poll() is None:
        try:
            tree = [root] + root.children(recursive=True)
        except psutil.NoSuchProcess:
            break
        rss = 0
        for p in tree:
            try:
                rss += p.memory_info().rss
                t = p.cpu_times()
                cpu_by_pid[p.pid] = t.user + t.system
            except psutil.Error:
                pass
        peak = max(peak, rss)
        time.sleep(0.1)
    stderr = proc.stderr.read().decode()
    wall = time.time() - start
    if proc.returncode != 0:
        sys.exit(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{stderr}")
    return wall, peak / 2**20, 100 * sum(cpu_by_pid.values()) / wall


def benchmark_target(target, experiments, num_runs):
    rows = []
    py = [sys.executable, "-m", "tamipami.cli"]
    for exp in experiments:
        cont_reads = count_fastq_pairs(exp["cont1"])
        exp_reads = count_fastq_pairs(exp["exp1"])
        print(f"[{target.label}] experiment {exp['experiment']} ({exp['enzyme']})", flush=True)
        for i in range(1, num_runs + 1):
            workdir = Path(tempfile.mkdtemp(prefix="tamipami_bench_run_"))
            h5 = workdir / "out.h5"
            try:
                pt, pm, pc = run_measured(
                    py + ["process", "--cont1", exp["cont1"], "--cont2", exp["cont2"],
                          "--exp1", exp["exp1"], "--exp2", exp["exp2"],
                          "--library", exp["library"], "--outfile", str(h5)], target.env)
                qt, qm, qc = run_measured(
                    py + ["predict", "--input", str(h5), "--cutoff", cutoff_for(h5),
                          "--predict_out", str(workdir / "predict")], target.env)
            finally:
                shutil.rmtree(workdir, ignore_errors=True)
            rows.append({"Experiment": exp["experiment"], "Enzyme": exp["enzyme"], "Run": i,
                         "process_time": pt, "predict_time": qt,
                         "process_max_mem": pm, "predict_max_mem": qm,
                         "process_avg_cpu": pc, "predict_avg_cpu": qc,
                         "Control_reads": cont_reads, "Exp_reads": exp_reads})
            print(f"  run {i}/{num_runs}: process {pt:.1f}s/{pm:.0f}MB  predict {qt:.1f}s/{qm:.0f}MB", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--num_runs", type=int, default=10, help="runs per experiment (default 10)")
    ap.add_argument("--experiments", default=str(bc.DEFAULT_EXPERIMENTS))
    ap.add_argument("--data-dir", default=str(bc.DEFAULT_DATA_DIR))
    ap.add_argument("--ref", action="append",
                    help=f"git tag/branch/commit or {bc.WORKTREE} (repeatable; default: latest tag and {bc.WORKTREE})")
    ap.add_argument("--outdir", default=str(bc.RESULTS_DIR))
    args = ap.parse_args()

    experiments = bc.read_experiments(args.experiments, args.data_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for ref in args.ref or bc.default_refs():
        with bc.resolve_target(ref) as target:
            if target.dirty:
                print(f"NOTE: {ref} has uncommitted changes; results are labeled {target.label}")
            rows = benchmark_target(target, experiments, args.num_runs)
            out = outdir / f"cli_{target.label}.csv"
            with open(out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS)
                w.writeheader()
                w.writerows(rows)
            bc.write_meta(out.with_suffix(".meta.json"), target,
                          {"num_runs": args.num_runs, "cutoff_value_all_lengths": CUTOFF_VALUE, "experiments": args.experiments})
            print(f"Wrote {out}")


if __name__ == "__main__":
    main()
