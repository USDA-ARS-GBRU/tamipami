# Degenerate-solver benchmarks

These benchmarks test `tamipami.degenerate.seqs_to_degenerates`, which finds the
smallest exact set of IUPAC degenerate sequences covering a selected set of PAM
k-mers. They answer two questions:

1. **Is the result globally optimal?** `run_oracle_benchmark.py` compares TAMIPAMI
   with an independent oracle (`oracle.py`). The oracle enumerates all 15^k IUPAC
   codes, keeps every code whose expansion lies inside the selected set, and solves
   the same three-level lexicographic exact-cover problem (patterns, then degenerate
   positions, then total bases) to proven optimality. It shares no code with
   TAMIPAMI's candidate generation.
2. **How do time, memory and quality change with the size of the selection?**
   `run_selection_sweep.py` ranks the k-mers of a real experiment by z-score and
   solves the top 5–100% at k = 3, 4 and 5. It can run against an older checkout to
   compare versions.

## Files

| File | Purpose |
|---|---|
| `oracle.py` | Independent oracle, plus `tamipami_solve`, which runs TAMIPAMI, checks that the result is an exact cover, and records whether every CP-SAT stage was proven optimal and whether candidate enumeration was capped |
| `run_oracle_benchmark.py` | Runs TAMIPAMI and the oracle on synthetic, random or real selections and appends one JSON line per set |
| `run_selection_sweep.py` | Times each selection size in its own process (peak memory per solve, 300 s kill timeout) |
| `summarize.py` | Writes `results/oracle_summary.csv` and `results/selection_summary.csv` from the JSON-lines files |
| `results/example_kmers.h5` | `tamipami process` output for the bundled SpCas9 example (`tamipami/assets/Example.zip`) |
| `results/oracle_*.jsonl` | Per-set oracle comparisons (2,365 sets) |
| `results/selection_current.jsonl` | Selection sweep with the current code (v0.4.0) |
| `results/selection_101284e.jsonl` | The same sweep with the code before the release audit (`101284e`) |
| `results/*_summary.csv` | Summary tables produced by `summarize.py` |

### Result fields

Scores are `[patterns, degenerate positions, bases]` and are compared
lexicographically (lower is better).

- **Oracle files:** `tamipami` and `oracle` are the two scores. `oracle` is `null`
  when the oracle could not prove optimality within 300 s per level. `match` means
  the scores are identical. `tamipami_proven` means all CP-SAT stages returned
  OPTIMAL and the candidate pool was not capped. `oracle_exact_codes` is the size
  of the oracle's full candidate pool.
- **Selection files:** `stage_statuses` lists the CP-SAT status of each stage.
  There are three stages, or six when the tree-only solve runs before the full
  solve. `peak_rss_mb` is the peak resident memory of the worker process. `codes`
  shows at most 10 codes.

## Reproducing

Run the commands from this directory, in the `tamipami` conda environment with the
package installed.

```bash
# 1. Example k-mer table (results/example_kmers.h5)
mkdir -p /tmp/tp_example && unzip -o ../../tamipami/assets/Example.zip -d /tmp/tp_example
tamipami process -c /tmp/tp_example/example_cont1.fastq.gz -c2 /tmp/tp_example/example_cont2.fastq.gz \
    -e /tmp/tp_example/example_exp1.fastq.gz -e2 /tmp/tp_example/example_exp2.fastq.gz \
    --library RTW554 --length 5 --outfile results/example_kmers.h5

# 2. Oracle comparison (about 2 h on 8 cores; dense random k=4 sets dominate)
python run_oracle_benchmark.py --example results/example_kmers.h5 --group "SpCas9 example" \
    --out results/oracle_example.jsonl
python run_oracle_benchmark.py --instances ../degenerate_multi/bench_multi_out/k3_m8/instances.jsonl \
    --group "synthetic k=3, 8 patterns" --out results/oracle_synthetic_k3_m8.jsonl
python run_oracle_benchmark.py --instances ../degenerate_multi/bench_k5_m4/k5_m4/instances.jsonl --limit 100 \
    --group "synthetic k=5, 4 patterns" --out results/oracle_synthetic_k5_m4.jsonl
python run_oracle_benchmark.py --instances ../degenerate_multi/bench_k4_m4/k3_m4/instances.jsonl \
    --group "synthetic k=3, 4 patterns" --out results/oracle_synthetic_k3_m4.jsonl
python run_oracle_benchmark.py --instances ../degenerate_multi/bench_k4_m4/k4_m4/instances.jsonl \
    --group "synthetic k=4, 4 patterns" --out results/oracle_synthetic_k4_m4.jsonl
python run_oracle_benchmark.py --random 3 --per-fraction 20 --out results/oracle_random_k3.jsonl
python run_oracle_benchmark.py --random 4 --per-fraction 5 --out results/oracle_random_k4.jsonl

# 3. Selection sweep, current code and the pre-audit commit
python run_selection_sweep.py --hdf results/example_kmers.h5 --label "v0.4.0 (current)" \
    --out results/selection_current.jsonl
git worktree add /tmp/tp_101284e 101284e
cp ../../tamipami/_version.py /tmp/tp_101284e/tamipami/   # generated file, not in git
python run_selection_sweep.py --hdf results/example_kmers.h5 --pythonpath /tmp/tp_101284e \
    --label "101284e (before audit)" --out results/selection_101284e.jsonl
git worktree remove /tmp/tp_101284e

# 4. Summary tables
python summarize.py results
```

`run_oracle_benchmark.py` appends to its output, so delete old `oracle_*.jsonl`
files before a full rerun. `run_selection_sweep.py` overwrites its output. The solver
uses CP-SAT's deterministic time limit, so repeated runs on the same machine give the
same covers. Wall-clock times depend on the machine. These results came from an
Apple-silicon Mac with `TAMIPAMI_SOLVER_WORKERS` unset (8 workers).

## Results

### Which selections are realistic

Even the most PAM-flexible engineered nucleases do not recognize every k-mer. SpRY,
for example, still prefers NRN over NYN, which already excludes about half of all
sequences. In practice, the PAM set passed to degenerate analysis is the automatic
cutoff or, at most, the top 25–50% of k-mers by z-score. Selections above 50%
appear below only to find the solver's limits.

### Global optimality (oracle comparison)

Of the 2,365 sets, the oracle proved the global optimum for 2,354. **TAMIPAMI's
cover matched the proven optimum at all three objective levels in all 2,354 sets.**
In 2,353 of them, TAMIPAMI also proved optimality itself. In the remaining set (a
random 90% selection of 3-mers), TAMIPAMI reached its budget but had already found
the optimum. The candidate cap was never reached; the largest candidate pool had
13,641 codes.

| Test sets | Sets | Sequences per set | Oracle proved optimum | TAMIPAMI = optimum |
|---|---|---|---|---|
| Synthetic ground truth, k=3, 4 patterns | 1,000 | 5–22 | 1,000 | 1,000 |
| Synthetic ground truth, k=3, 8 patterns | 100 | 13–29 | 100 | 100 |
| Synthetic ground truth, k=4, 4 patterns | 1,000 | 5–19 | 1,000 | 1,000 |
| Synthetic ground truth, k=5, 4 patterns | 100 | 4–11 | 100 | 100 |
| Random k=3 subsets, 5–90% | 120 | 3–58 | 120 | 120 |
| Random k=4 subsets, 5–90% | 30 | 13–230 | 20 | 20 |
| SpCas9 example, auto cutoff and top 5–50%, k=3–5 | 15 | 3–512 | 14 | 14 |
| **Total** | **2,365** | | **2,354** | **2,354** |

The oracle did not finish on 11 sets. All of them are dense: the SpCas9 top 50% of
5-mers (512 sequences) and the random k=4 sets at 75% and 90%. TAMIPAMI returned a
valid exact cover for each of them in under 7 s, but without a proof of minimality.

**By selection size, for realistic selections (≤ 50%):**

- **Automatic cutoff and ≤ 25%** (all lengths): every set was solved to proven
  global optimality. The SpCas9 auto-cutoff selections (NGG, NGGN, NGGNN; 4–64
  sequences) take under 0.05 s, and the top 25% of 5-mers (256 sequences) takes 7 s.
- **50%:** every k=3 and k=4 set was solved to proven optimality (k=4, 128
  sequences: 4 s). The single exception is the SpCas9 top 50% of 5-mers. That set
  has no proven optimum from either TAMIPAMI or the oracle, and TAMIPAMI returns
  a feasible cover in 6 s.

### Time, memory and quality by selection size (SpCas9 example)

`results/selection_summary.csv`; scores are patterns/degenerate positions/bases.

| k | Selection | n | Before audit (101284e) | Current | Current proven | Current time (s) | Current peak RSS (MB) |
|---|---|---|---|---|---|---|---|
| 3 | auto | 4 | 1/1/6 | 1/1/6 (NGG) | yes | 0.02 | 279 |
| 3 | 25% | 16 | 8/5/31 | 8/5/31 | yes | 0.01 | 289 |
| 3 | 50% | 32 | 11/12/48 | 11/12/48 | yes | 0.4 | 304 |
| 4 | auto | 16 | 1/2/10 | 1/2/10 (NGGN) | yes | 0.01 | 291 |
| 4 | 25% | 64 | 23/20/121 | 23/20/121 | yes | 0.6 | 326 |
| 4 | 50% | 128 | 35/47/202 | 35/47/202 | yes | 4.3 | 449 |
| 4 | 75% | 192 | 33/60/214 | 33/60/216 | no | 4.5 | 730 |
| 4 | 90% | 230 | 24/47/168 | 26/45/177 | no | 5.1 | 853 |
| 4 | 100% | 256 | 1/4/16 (unproven) | 1/4/16 (NNNN) | yes | 0.04 | 278 |
| 5 | auto | 64 | 1/3/14 | 1/3/14 (NGGNN) | yes | 0.01 | 267 |
| 5 | 10% | 102 | 32/10/176 | 32/10/176 | yes | 2.5 | 361 |
| 5 | 25% | 256 | 102/80/604 | 102/80/604 | yes | 7.7 | 621 |
| 5 | 50% | 512 | 128/187/863 | 130/191/874 | no | 6.5 | 810 |
| 5 | 75% | 768 | 150/199/1005 | 168/218/1206 | no | 8.5 | 1,073 |
| 5 | 90% | 922 | 922 singletons | 270/280/1796 | no (pool capped) | 8.9 | 1,106 |
| 5 | 100% | 1024 | crashed (RuntimeError) | 1/5/20 (NNNNN) | yes | 0.3 | 302 |

- Up to 25% at every length, plus 50% at k ≤ 4, the two versions give identical,
  proven-optimal covers. The current code is sometimes slower on small sets because
  it solves the tree-only candidates first, but no set in this range takes more
  than 8 s.
- Above 50% at k ≥ 4, neither version proves optimality. The pre-audit code used a
  wall-clock limit with every core, so its results depended on the machine and
  sometimes beat the current code's (k=5 75%: 150 vs 168 patterns). The current code
  uses a fixed deterministic budget, so its results are reproducible. It never
  returns 922 singletons (the old k=5 90% result) and does not crash on a full
  selection (the old k=5 100% result).
- Peak memory stays under 1 GB for realistic selections and reaches 1.1 GB at the
  densest k=5 selections.

**Practical implication.** Restricting set analysis to the automatic cutoff, or to
at most the top 25–50% of k-mers, keeps every tested selection except one within
the range where TAMIPAMI proves a globally optimal cover in seconds. That exception,
the top 50% of 5-mers, is a selection that a real nuclease is unlikely to produce.
