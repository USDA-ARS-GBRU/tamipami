import time
import csv
import psutil
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Configuration
STREAMLIT_URL = "http://localhost:8501"
N_RUNS = 10
# DATA_FILES = {
#     "cont1": "/Users/adam.rivers/Documents/tamipami/data/spcas9_RTW554/SpCas9_cont_R1_001.fastq.gz",
#     "cont2": "/Users/adam.rivers/Documents/tamipami/data/spcas9_RTW554/SpCas9_cont_R2_001.fastq.gz",
#     "exp1": "/Users/adam.rivers/Documents/tamipami/data/spcas9_RTW554/SpCas9_exp_R1_001.fastq.gz",
#     "exp2": "/Users/adam.rivers/Documents/tamipami//data/spcas9_RTW554/SpCas9_exp_R2_001.fastq.gz"
# }
# LIBRARY = "RTW554"
# CSV_FILE = "/Users/adam.rivers/Documents/tamipami/benchmark/benchmark_app_chrome_spcas9_apple_m2pro.csv"
DATA_FILES = {
    "cont1": "/Users/adam.rivers/Documents/tamipami/data/cas12i_RTW574/Cas12i_cont_R1_001.fastq.gz",
    "cont2": "/Users/adam.rivers/Documents/tamipami/data/cas12i_RTW574/Cas12i_cont_R2_001.fastq.gz",
    "exp1": "/Users/adam.rivers/Documents/tamipami/data/cas12i_RTW574/Cas12i1_exp_R1_001.fastq.gz",
    "exp2": "/Users/adam.rivers/Documents/tamipami/data/cas12i_RTW574/Cas12i1_exp_R2_001.fastq.gz"
}
LIBRARY = "RTW574"
CSV_FILE = "/Users/adam.rivers/Documents/tamipami/benchmark/benchmark_app_chrome_cas12i_apple_m2pro.csv"
def run_once(driver, run_number):
    driver.get(STREAMLIT_URL)
    wait = WebDriverWait(driver, 60)

    # Upload files to the correct uploaders
    uploaders = wait.until(EC.presence_of_all_elements_located((By.XPATH, "//input[@type='file']")))
    for idx, (key, path) in enumerate(DATA_FILES.items()):
        uploaders[idx].send_keys(path)
        time.sleep(1)

    # Open the selectbox
    selectbox_div = wait.until(EC.element_to_be_clickable((By.XPATH, "//div[contains(@data-baseweb, 'select')]")))
    selectbox_div.click()
    time.sleep(0.5)

    # Select the option "RTW554"
    option = wait.until(EC.element_to_be_clickable(
        (By.XPATH, "//div[contains(@class, 'st-emotion-cache-qiev7j') and text()='RTW554']")))
    option.click()
    time.sleep(0.5)

    # Wait for the submit button to be enabled and click it
    # Wait for the submit button to be clickable and click it
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
    import sys
    chrome_options = Options()
    chrome_options.add_argument("--start-maximized")
    driver = webdriver.Chrome(options=chrome_options)

    results = []
    for i in range(N_RUNS):
        print(f"Run {i+1}/{N_RUNS}...")
        elapsed, max_mem, avg_cpu = run_once(driver, i+1)
        results.append({"run": i+1, "elapsed_time_sec": elapsed, "max_mem_mb": max_mem, "avg_cpu_percent": avg_cpu})
        print(f"  Time: {elapsed:.2f}s, Max Mem: {max_mem:.1f}MB, Avg CPU: {avg_cpu:.1f}%")
        time.sleep(2)  # Optional: wait between runs

    driver.quit()

    # Write results to CSV
    with open(CSV_FILE, "w", newline="") as csvfile:
        fieldnames = ["run", "elapsed_time_sec", "max_mem_mb", "avg_cpu_percent"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nResults written to {CSV_FILE}")

if __name__ == "__main__":
    main()