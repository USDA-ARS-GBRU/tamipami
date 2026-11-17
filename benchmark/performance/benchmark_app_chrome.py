#!/usr/bin/env python3
"""Benchmark the TamiPami Streamlit app by driving it with Chrome (Selenium).

For each target ref (default: latest release tag and the working tree) this script starts its own
``tamipami serve`` on a free port, submits every experiment N times through the real UI, and records
wall time plus CPU/memory of the server process tree (where the work happens) and of the browser.
Results go to ``results/app_<label>.csv`` with a ``.meta.json`` sidecar.

Requires: selenium (with Chrome installed; Selenium Manager fetches the driver) and psutil.
"""
import argparse
import csv
import os
import socket
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path

import psutil
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

import bench_common as bc

FIELDS = ["experiment", "enzyme", "library", "run", "elapsed_time_sec", "server_baseline_mem_mb", "server_max_mem_mb",
          "server_avg_cpu_percent", "browser_max_mem_mb", "browser_avg_cpu_percent"]
RESULT_TIMEOUT_S = 1800


class TreeSampler(threading.Thread):
    """Samples RSS and cumulative CPU time of a process and all its descendants."""

    def __init__(self, root_pid):
        super().__init__(daemon=True)
        self.root = psutil.Process(root_pid)
        self.stop = threading.Event()
        self.peak = 0
        self.cpu_by_pid = {}
        self.cpu_base = {}
        self.primed = False
        self.baseline_rss = 0
        for p in [self.root] + self.root.children(recursive=True):
            try:
                self.baseline_rss += p.memory_info().rss
            except psutil.Error:  # exited or zombie (Chrome leaves these behind)
                pass

    def run(self):
        while not self.stop.is_set():
            try:
                tree = [self.root] + self.root.children(recursive=True)
            except psutil.Error:
                break
            rss = 0
            for p in tree:
                try:
                    rss += p.memory_info().rss
                    t = p.cpu_times()
                    cpu = t.user + t.system
                    # cpu_times is cumulative since process start; keep only the increase while sampling
                    self.cpu_base.setdefault(p.pid, 0.0 if self.primed else cpu)
                    self.cpu_by_pid[p.pid] = cpu - self.cpu_base[p.pid]
                except psutil.Error:
                    pass
            self.peak = max(self.peak, rss)
            self.primed = True
            time.sleep(0.1)

    def finish(self, wall):
        self.stop.set()
        self.join()
        return self.peak / 2**20, 100 * sum(self.cpu_by_pid.values()) / wall


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def start_server(target):
    port = free_port()
    env = dict(target.env, STREAMLIT_SERVER_HEADLESS="true",
               STREAMLIT_BROWSER_GATHER_USAGE_STATS="false")
    proc = subprocess.Popen(
        [sys.executable, "-m", "tamipami.cli", "serve", "--port", str(port), "--host", "localhost"],
        env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL,
        start_new_session=True)
    url = f"http://localhost:{port}"
    for _ in range(120):
        try:
            if urllib.request.urlopen(url + "/_stcore/health", timeout=2).status == 200:
                return proc, url
        except Exception:
            time.sleep(0.5)
    proc.kill()
    raise RuntimeError("Streamlit server did not become healthy")


def stop_server(proc):
    try:
        for child in psutil.Process(proc.pid).children(recursive=True):
            child.kill()
    except psutil.Error:
        pass
    proc.kill()
    proc.wait()


def run_once(driver, url, server_pid, exp):
    driver.get(url)
    wait = WebDriverWait(driver, 60)
    uploaders = wait.until(EC.presence_of_all_elements_located((By.XPATH, "//input[@type='file']")))
    for idx, key in enumerate(("cont1", "cont2", "exp1", "exp2")):
        uploaders[idx].send_keys(str(Path(exp[key]).resolve()))
        time.sleep(1)

    wait.until(EC.element_to_be_clickable((By.XPATH, "//div[contains(@data-baseweb, 'select')]"))).click()
    time.sleep(0.5)
    wait.until(EC.element_to_be_clickable(
        (By.XPATH, f"//*[@role='option'][.//*[normalize-space(text())='{exp['library']}'] "
                   f"or normalize-space(text())='{exp['library']}']"))).click()
    time.sleep(0.5)
    submit = wait.until(EC.element_to_be_clickable(
        (By.XPATH, "//button[@data-testid='stBaseButton-secondaryFormSubmit']")))

    server = TreeSampler(server_pid)
    browser = TreeSampler(driver.service.process.pid)  # chromedriver and the Chrome it spawned
    server.start()
    browser.start()
    start = time.time()
    submit.click()
    done = False
    while time.time() - start < RESULT_TIMEOUT_S:
        if driver.find_elements(By.XPATH, "//*[contains(text(), 'Library information')]"):
            done = True
            break
        time.sleep(0.1)
    elapsed = time.time() - start
    s_mem, s_cpu = server.finish(elapsed)
    b_mem, b_cpu = browser.finish(elapsed)
    if not done:
        raise RuntimeError(f"Timed out waiting for results after {RESULT_TIMEOUT_S}s")
    return elapsed, server.baseline_rss / 2**20, s_mem, s_cpu, b_mem, b_cpu


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--num-runs", "--num_runs", dest="num_runs", type=int, default=10)
    ap.add_argument("--experiments", default=str(bc.DEFAULT_EXPERIMENTS))
    ap.add_argument("--data-dir", default=str(bc.DEFAULT_DATA_DIR))
    ap.add_argument("--ref", action="append",
                    help=f"git tag/branch/commit or {bc.WORKTREE} (repeatable; default: latest tag and {bc.WORKTREE})")
    ap.add_argument("--outdir", default=str(bc.RESULTS_DIR))
    ap.add_argument("--headless", action="store_true", help="run Chrome headless")
    args = ap.parse_args()

    experiments = bc.read_experiments(args.experiments, args.data_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for ref in args.ref or bc.default_refs():
        with bc.resolve_target(ref) as target:
            opts = Options()
            opts.add_argument("--headless=new" if args.headless else "--start-maximized")
            driver = webdriver.Chrome(options=opts)
            rows = []
            server = None
            try:
                for exp in experiments:
                    print(f"[{target.label}] experiment {exp['experiment']} ({exp['enzyme']})", flush=True)
                    for i in range(1, args.num_runs + 1):
                        # Fresh server per run: a long-lived one keeps memory from earlier sessions
                        # (baseline grew 0.8 -> 3.6 GB over 7 runs), which would contaminate the peaks.
                        server, url = start_server(target)
                        time.sleep(2)
                        t, s0, sm, sc, bm, bcpu = run_once(driver, url, server.pid, exp)
                        rows.append({"experiment": exp["experiment"], "enzyme": exp["enzyme"],
                                     "library": exp["library"], "run": i, "elapsed_time_sec": t,
                                     "server_baseline_mem_mb": s0, "server_max_mem_mb": sm, "server_avg_cpu_percent": sc,
                                     "browser_max_mem_mb": bm, "browser_avg_cpu_percent": bcpu})
                        print(f"  run {i}/{args.num_runs}: {t:.1f}s server {sm:.0f}MB", flush=True)
                        stop_server(server)
                        server = None
                        time.sleep(2)
            finally:
                driver.quit()
                if server:
                    stop_server(server)
            out = outdir / f"app_{target.label}.csv"
            with open(out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS)
                w.writeheader()
                w.writerows(rows)
            bc.write_meta(out.with_suffix(".meta.json"), target,
                          {"num_runs": args.num_runs, "headless": args.headless, "cpu_delta": True,
                           "experiments": args.experiments})
            print(f"Wrote {out}")


if __name__ == "__main__":
    main()
