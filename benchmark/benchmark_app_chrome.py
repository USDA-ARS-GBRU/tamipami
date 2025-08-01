import time
import csv
import psutil
import argparse
import os
import sys
import tempfile
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

STREAMLIT_URL = "http://localhost:8501"

def run_once(driver, data_files, library):
    driver.get(STREAMLIT_URL)
    wait = WebDriverWait(driver, 60)

    # Upload files to the correct uploaders
    uploaders = wait.until(EC.presence_of_all_elements_located((By.XPATH, "//input[@type='file']")))
    for idx, (key, path) in enumerate(data_files.items()):
        uploaders[idx].send_keys(path)
        time.sleep(1)

    # Open the selectbox
    selectbox_div = wait.until(EC.element_to_be_clickable((By.XPATH, "//div[contains(@data-baseweb, 'select')]")))
    selectbox_div.click()
    time.sleep(0.5)

    # Select the option for the library
    option = wait.until(EC.element_to_be_clickable(
        (By.XPATH, f"//div[contains(@class, 'st-emotion-cache-qiev7j') and text()='{library}']")))
    option.click()
    time.sleep(0.5)

    # Wait for the submit button to be enabled and click it
    submit_btn = wait.until(EC.element_to_be_clickable((By.XPATH, "//button[@data-testid='stBaseButton-secondaryFormSubmit']")))
    chrome_pid = driver.service.process.pid
    chrome_proc = psutil.Process(chrome_pid)

    # Start monitoring
    cpu_samples = []
    mem_samples = []
    start = time.time()
    submit_btn.click()

    # Wait for results
    while True:
        try:
            if driver.find_elements(By.XPATH, "//*[contains(text(), 'Library information')]"):
                break
            cpu_samples.append(chrome_proc.cpu_percent(interval=0.1))
            mem_samples.append(chrome_proc.memory_info().rss)
        except Exception:
            break

    elapsed = time.time() - start
    max_mem = max(mem_samples) / (1024 * 1024) if mem_samples else 0  # MB
    avg_cpu = sum(cpu_samples) / len(cpu_samples) if cpu_samples else 0  # %

    return elapsed, max_mem, avg_cpu

def main():
    parser = argparse.ArgumentParser(description="Benchmark Streamlit app with Chrome automation.")
    parser.add_argument("--experiments", type=str, required=True, help="CSV file with experiment definitions")
    parser.add_argument("--data-dir", type=str, default="", help="Directory where experiment data files are located")
    parser.add_argument("-n", "--num_runs", type=int, default=3, help="Number of runs per experiment")
    parser.add_argument("--output", type=str, default="benchmark_app_chrome_results.csv", help="Output CSV file")
    args = parser.parse_args()

    # Read experiments CSV
    with open(args.experiments, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        experiments = list(reader)

    def normkey(k):
        return k.strip().lower().replace('\ufeff', '')
    if experiments:
        keymap = {normkey(k): k for k in experiments[0].keys()}
    else:
        keymap = {}

    chrome_options = Options()
    chrome_options.add_argument("--start-maximized")
    driver = webdriver.Chrome(options=chrome_options)

    results = []
    for exp in experiments:
        exp_num = exp.get(keymap.get('experiment', 'Experiment'), '')
        enzyme = exp.get(keymap.get('enzyme', 'Enzyme'), '')
        library = exp.get(keymap.get('control library', 'Control Library'), '')
        data_dir = args.data_dir
        def join_data_dir(filename):
            return os.path.join(data_dir, filename) if data_dir and filename else filename
        data_files = {
            "cont1": join_data_dir(exp.get(keymap.get('cont1', 'cont1'), '')),
            "cont2": join_data_dir(exp.get(keymap.get('cont2', 'cont2'), '')),
            "exp1": join_data_dir(exp.get(keymap.get('exp1', 'exp1'), '')),
            "exp2": join_data_dir(exp.get(keymap.get('exp2', 'exp2'), '')),
        }
        print(f"Running experiment {exp_num} ({enzyme})...")
        for i in range(args.num_runs):
            print(f"Run {i+1}/{args.num_runs}...")
            elapsed, max_mem, avg_cpu = run_once(driver, data_files, library)
            results.append({
                "experiment": exp_num,
                "enzyme": enzyme,
                "library": library,
                "run": i+1,
                "elapsed_time_sec": elapsed,
                "max_mem_mb": max_mem,
                "avg_cpu_percent": avg_cpu
            })
            print(f"  Time: {elapsed:.2f}s, Max Mem: {max_mem:.1f}MB, Avg CPU: {avg_cpu:.1f}%")
            time.sleep(2)  # Optional: wait between runs

    driver.quit()

    # Write results to CSV
    with open(args.output, "w", newline="") as csvfile:
        fieldnames = ["experiment", "enzyme", "library", "run", "elapsed_time_sec", "max_mem_mb", "avg_cpu_percent"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nResults written to {args.output}")

if __name__ == "__main__":
    main()