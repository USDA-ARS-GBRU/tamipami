#!/usr/bin/env python3
"""
Benchmark script for tamipami CLI: runs `process` and `predict` N times, and reports runtime, CPU, and memory statistics.
"""
import argparse
import subprocess
import time
import psutil
import shutil
import os
import sys
import tempfile
from statistics import mean
import csv
import gzip

CUTOFF = '{"3":2,"4":2,"5":2,"6":2}'

def count_gzipped_fastq_reads(filepath):
    """
    Reads a gzipped file, uncompresses it, counts the lines,
    and returns the line count divided by 4.

    Args:
        filepath (str): The path to the gzipped file.

    Returns:
        float: The number of lines in the uncompressed file divided by 4.
               Returns 0 if the file cannot be read or is empty.
    """
    line_count = 0
    try:
        with gzip.open(filepath, 'rt') as f:
            for _ in f:
                line_count += 1
        return line_count / 4
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return 0
    except Exception as e:
        print(f"An error occurred: {e}")
        return 0
    

def run_once(process_cmd, predict_cmd):
    # Run process, write output to file, measure resource usage
    start_time = time.time()
    proc = subprocess.Popen(process_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps_proc = psutil.Process(proc.pid)
    max_mem = 0
    cpu_times = []
    while proc.poll() is None:
        try:
            mem = ps_proc.memory_info().rss
            max_mem = max(max_mem, mem)
            cpu_times.append(ps_proc.cpu_percent(interval=0.1))
        except psutil.NoSuchProcess:
            break
    _, process_stderr = proc.communicate()
    process_time = time.time() - start_time
    if proc.returncode != 0:
        print(f"\nError running process command: {' '.join(process_cmd)}", file=sys.stderr)
        print(process_stderr.decode(), file=sys.stderr)
        sys.exit(proc.returncode)

    # Run predict, read from file, measure resource usage
    start_time2 = time.time()
    proc2 = subprocess.Popen(predict_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps_proc2 = psutil.Process(proc2.pid)
    max_mem2 = 0
    cpu_times2 = []
    while proc2.poll() is None:
        try:
            mem2 = ps_proc2.memory_info().rss
            max_mem2 = max(max_mem2, mem2)
            cpu_times2.append(ps_proc2.cpu_percent(interval=0.1))
        except psutil.NoSuchProcess:
            break
    _, predict_stderr = proc2.communicate()
    predict_time = time.time() - start_time2
    if proc2.returncode != 0:
        print(f"\nError running predict command: {' '.join(predict_cmd)}", file=sys.stderr)
        print(predict_stderr.decode(), file=sys.stderr)
        sys.exit(proc2.returncode)

    return {
        "process_time": process_time,
        "predict_time": predict_time,
        "process_max_mem": max_mem / (1024 * 1024),
        "predict_max_mem": max_mem2 / (1024 * 1024),
        "process_avg_cpu": mean(cpu_times) if cpu_times else 0,
        "predict_avg_cpu": mean(cpu_times2) if cpu_times2 else 0,
    }

def main():
    parser = argparse.ArgumentParser(description="Benchmark tamipami CLI process and predict commands.")
    parser.add_argument("-n", "--num_runs", type=int, default=3, help="Number of runs (N)")
    parser.add_argument("--experiments", type=str, help="CSV file with experiment definitions")
    parser.add_argument("--data-dir", type=str, default="", help="Directory where experiment data files are located")
    args = parser.parse_args()

    results = []
    if args.experiments:
        # Read experiments CSV
        with open(args.experiments, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            experiments = list(reader)

        # Normalize header keys for robust access
        def normkey(k):
            return k.strip().lower().replace('\ufeff', '')
        # Build a mapping from normalized key to actual key
        if experiments:
            keymap = {normkey(k): k for k in experiments[0].keys()}
        else:
            keymap = {}
        for exp in experiments:
            # Use normalized keys for access
            exp_num = exp.get(keymap.get('experiment', 'Experiment'), '')
            enzyme = exp.get(keymap.get('enzyme', 'Enzyme'), '')
            data_dir = args.data_dir
            def join_data_dir(filename):
                return os.path.join(data_dir, filename) if data_dir and filename else filename
            cont1 = join_data_dir(exp.get(keymap.get('cont1', 'cont1'), ''))
            cont2 = join_data_dir(exp.get(keymap.get('cont2', 'cont2'), ''))
            exp1 = join_data_dir(exp.get(keymap.get('exp1', 'exp1'), ''))
            exp2 = join_data_dir(exp.get(keymap.get('exp2', 'exp2'), ''))
            control_lib = exp.get(keymap.get('control library', 'Control Library'), '')
            cont_read_cnt = count_gzipped_fastq_reads(cont1)
            exp_read_cnt = count_gzipped_fastq_reads(exp1)
            print(f"Running experiment {exp_num} ({enzyme})...")

            for i in range(args.num_runs):
                with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as process_outfile:
                    predict_outdir = tempfile.mkdtemp()
                    process_cmd = [sys.executable, "-m", "tamipami.cli", "process",
                        "--cont1", cont1, "--cont2", cont2, "--exp1", exp1, "--exp2", exp2,
                        "--library", control_lib, "--outfile", process_outfile.name]
                    predict_cmd = [sys.executable, "-m", "tamipami.cli", "predict",
                        "--input", process_outfile.name,
                        "--cutoff", CUTOFF, "--predict_out", predict_outdir]
                    print(f"Run {i+1}/{args.num_runs}...")
                    stats = run_once(process_cmd, predict_cmd)
                    stats["Experiment"] = exp_num
                    stats["Enzyme"] = enzyme
                    stats["Control_reads"] = cont_read_cnt
                    stats["Exp_reads"] = exp_read_cnt
                    results.append(stats)
                    print(f"  process_time: {stats['process_time']:.2f}s, predict_time: {stats['predict_time']:.2f}s, "
                        f"process_max_mem: {stats['process_max_mem']:.1f}MB, predict_max_mem: {stats['predict_max_mem']:.1f}MB, "
                        f"process_avg_cpu: {stats['process_avg_cpu']:.1f}%, predict_avg_cpu: {stats['predict_avg_cpu']:.1f}%")
                    os.remove(process_outfile.name)
                    shutil.rmtree(predict_outdir)
    # Write results to CSV
    if results:
        csv_file = "benchmark_results_v014.csv"
        fieldnames = list(results[0].keys())
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f"\nResults written to {csv_file}")

        # Print summary
        print("\nSummary (averages over runs):")
        keys = [k for k in results[0].keys() if k not in ("Experiment", "Enzyme")]
        for k in keys:
            avg = mean([r[k] for r in results])
            print(f"{k}: {avg:.2f}")

if __name__ == "__main__":
    main()
