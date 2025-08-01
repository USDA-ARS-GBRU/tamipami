# Benchmarking performance of Tamipami

## command line benchmarking

The script `run_benchmark.py` runs each experimental dataset in the paper N times and reports the timing and resource use in an output CSV.

```{bash}
python run_benchmark.py --experiments ../data/metadata/experiments.csv --data-dir ../data/fastq --num_runs 10
```
## Web application benchmarking

The script `benchmark_app_chrome.py` runs experimental dataset in the paper through a locally hosted instance of the Tamipami web app N times and reports the timing and resource use in an output CSV.  It uses the Python  browser testing framework [Selenium](https://www.selenium.dev/documentation/) to automatically interact with the web application and measure client-side resource use and timings.

Streamlit, the web application framework for Tamipami, uses cashing to improve performance. To simulate the first run experience, we disable caching by commenting out the decorator `@st.cache_data` above the `process()` function of `app.py`.

```{bash}
python benchmark_app_chrome.py --experiments /Users/adam.rivers/Documents/tamipami/data/metadata/experiments.csv --data-dir /Users/adam.rivers/Documents/tamipami/data/fastq --num_runs 1 --output chrome_app_results.csv
```

## Test system configuration

Performance testing was conducted on a Macbook Pro with the following configuration:

Attribute | Value
----------|------
Model Name | MacBook Pro
Model Identifier | Mac14,9
Model Number | MPHJ3LL/A
Chip | Apple M2 Pro
Total Number of Cores | 12 (8 performance and 4 efficiency)
Memory | 16 GB
System Firmware Version | 11881.121.1
OS Loader Version | 11881.121.1
MacOS Version | 15.5 (24F74)
Chrome Version | 138.0.7204.184 (Official Build) (arm64)
Python version | 3.13.3


## Results