# Multi-pattern degenerate benchmark

`benchmark_multi.py` generates sets of k-mers from a known set of non-overlapping
degenerate IUPAC sequences (the "ground truth", GT), asks the solver for the minimal
exact cover, and compares the two.

**This is not an accuracy benchmark.** The planted GT is not guaranteed to be the
smallest cover, so a result that differs from it is not an error. Accuracy against a
proven global optimum is measured by the oracle comparison in
[`../degenerate_solver/`](../degenerate_solver/README.md). This benchmark's
`instances.jsonl` files are inputs to that comparison, so keep them.

## What each instance does

1. Samples k-mers and tiles them into `m` GT rectangles by cutting a Hamming-distance
   dendrogram. Only subtrees whose expansion lies inside the sampled set are kept.
2. Takes the union of the GT expansions as the input set `D`.
3. Solves with `degenerate.solve_degenerates`, the call the app and CLI use, and
   records whether the result was proven optimal (`proven_optimal`).
4. Enumerates the candidate pool separately and checks that every GT code is in it
   (`gt_missing_in_candidates`). Each GT code is exact, so a complete pool must contain
   all of them. The GT is **not** added to the pool.
5. Checks that the solver's codes are disjoint and their union equals `D`
   (`coverage_ok`), then compares them with the GT.

## Reading the results

`coverage_ok` is the correctness check: the answer is an exact cover. It is true for
every saved instance.

The verdicts compare the **number of codes** with the GT:

| Verdict | Meaning |
|---|---|
| `equivalent` | Same number of codes and the same expansions as the GT (`exact_block_matches`) |
| `better` | Fewer codes than the GT, so the planted set was not minimal |
| `valid-different` | Same number of codes as the GT but different codes. Not an error. Both are minimal in count, and the solver then minimizes degenerate positions and total bases, so its set is often more compact |
| `worse` | More codes than the GT. Never seen |
| `invalid` | Not an exact cover. Never seen |

`exact_block_matches` and `valid-different` are **not** error rates. A low
`exact_block_matches` mostly reflects an ambiguous GT, not a solver failure. The
informative columns are `worse` and `invalid`, which are 0 in every saved run, and
`better`, which shows the planted set was not minimal.

Saved runs (counted from the `instances.jsonl` files). In every instance the solver's
result was proven optimal, and every GT code was in the candidate pool.

| Directory | k | m | Instances | equivalent | valid-different | better | worse | invalid |
|---|---|---|---|---|---|---|---|---|
| `bench_multi_out/k3_m8` | 3 | 8 | 100 | 3 | 48 | 49 | 0 | 0 |
| `bench_k4_m4/k3_m4` | 3 | 4 | 1000 | 447 | 499 | 54 | 0 | 0 |
| `bench_k4_m4/k4_m4` | 4 | 4 | 1000 | 605 | 376 | 19 | 0 | 0 |
| `bench_k5_m4/k5_m4` | 5 | 4 | 1000 | 799 | 197 | 4 | 0 | 0 |

In the `valid-different` instances, the solver's secondary score (degenerate
positions, then total bases) is never larger than the GT's.

The pool check (`gt_missing_in_candidates` empty in all 3,100 instances) is a direct
test of candidate-generation completeness. It only tests GT codes, though. The
oracle comparison in `../degenerate_solver/` is the broader test.

## Limitations

- The GT is a sample of rectangles built by one procedure, not an arbitrary hard set,
  and it is not guaranteed minimal. Accuracy against proven optima is measured by the
  oracle benchmark, not here.
- Each run writes `summary.csv` into its own `<outdir>/k<k>_m<m>/` folder, so runs
  that share an `--outdir` no longer overwrite each other.
- The instances in this folder were regenerated with the current script. They differ
  from earlier versions in a small fraction of instances (for example 16 of 1,000 at
  k=3, 4 patterns), so the oracle results were recomputed from them.

## Running

```bash
python benchmark_multi.py --k 3 --m 8 --samples 100 --outdir bench_multi_out
python benchmark_multi.py --k 4 --m 4 --samples 1000 --outdir bench_k4_m4
```

Instances are written to `<outdir>/k<k>_m<m>/instances.jsonl`, one JSON record per
instance, with the GT codes, solver codes, verdict and compactness sums. The master
seed is fixed (12345), so reruns give the same instances.
