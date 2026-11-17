#!/usr/bin/env python
"""Compare tamipami's degenerate covers with an independent global-optimum oracle.

Test sets (choose one per run):
  --instances FILE   Synthetic ground-truth instances (instances.jsonl from
                     ../degenerate_multi); D is the union of the ground-truth codes.
  --random K         Random subsets of all K-mers at 5, 10, 25, 50, 75 and 90%.
  --example HDF      Selections from a `tamipami process` HDF5 file at lengths 3-5:
                     the automatic cutoff and the top 5, 10, 25 and 50% by z-score.

Each set is solved by tamipami.degenerate.seqs_to_degenerates and by the oracle
(oracle.py); one JSON line per set is appended to --out.

Example:
  python run_oracle_benchmark.py --random 3 --per-fraction 20 --out results/oracle_random_k3.jsonl
"""
import argparse
import itertools
import json
import random
import time
import warnings

from oracle import expand, oracle, tamipami_solve

FRACTIONS = (0.05, 0.1, 0.25, 0.5, 0.75, 0.9)


def instance_sets(path, limit):
    with open(path) as fh:
        for line in itertools.islice(fh, limit):
            rec = json.loads(line)
            yield None, sorted({s for c in rec["gt_codes"] for s in expand(c)})


def random_sets(k, per_fraction, seed):
    rng = random.Random(seed)
    allk = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
    for frac in FRACTIONS:
        for _ in range(per_fraction):
            yield frac, sorted(rng.sample(allk, max(2, round(frac * len(allk)))))


def example_sets(hdf):
    from tamipami import pam, tpio

    data = tpio.import_hdf(filename=hdf)
    obj = pam.pamSeqExp(multikmerdict=data)
    for k in (3, 4, 5):
        df = data[k].sort_values("zscore", ascending=False)
        auto = df[df.zscore >= obj.find_breakpoint(k)]
        yield ("auto", sorted(auto.kmers))
        for frac in (0.05, 0.1, 0.25, 0.5):
            yield frac, sorted(df.kmers.iloc[: round(frac * len(df))])


def main():
    warnings.filterwarnings("ignore")
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--instances", help="instances.jsonl with gt_codes")
    src.add_argument("--random", type=int, metavar="K", help="random subsets of all K-mers")
    src.add_argument("--example", metavar="HDF", help="tamipami process HDF5 output")
    p.add_argument("--limit", type=int, default=None, help="max instances to read (--instances)")
    p.add_argument("--per-fraction", type=int, default=20, help="random sets per fraction (--random)")
    p.add_argument("--seed", type=int, default=42, help="random seed (--random)")
    p.add_argument("--group", help="label for this group of sets (default: derived from the source)")
    p.add_argument("--oracle-time-limit", type=float, default=300.0, help="seconds per objective level")
    p.add_argument("--out", required=True, help="JSON-lines output file (appended)")
    args = p.parse_args()

    if args.instances:
        group = args.group or args.instances
        sets = instance_sets(args.instances, args.limit)
    elif args.random:
        group = args.group or f"random k={args.random}"
        sets = random_sets(args.random, args.per_fraction, args.seed)
    else:
        group = args.group or "example"
        sets = example_sets(args.example)

    with open(args.out, "a") as out:
        for fraction, D in sets:
            if len(D) < 2:
                continue
            t0 = time.time()
            tp_score, proven, limited = tamipami_solve(D)
            t1 = time.time()
            orc_score, n_exact = oracle(D, time_limit=args.oracle_time_limit)
            t2 = time.time()
            rec = dict(
                group=group, k=len(D[0]), fraction=fraction, n=len(D),
                tamipami=tp_score, tamipami_proven=proven, candidates_limited=limited,
                oracle=orc_score, oracle_exact_codes=n_exact,
                match=(orc_score is not None and tuple(tp_score) == tuple(orc_score)),
                t_tamipami=round(t1 - t0, 2), t_oracle=round(t2 - t1, 2),
            )
            out.write(json.dumps(rec) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
