#!/usr/bin/env python
"""Summarize oracle and selection-sweep results into CSV tables.

  python summarize.py results/

writes results/oracle_summary.csv and results/selection_summary.csv and prints both.
"""
import csv
import glob
import json
import os
import sys
from collections import OrderedDict


def load(pattern):
    rows = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as fh:
            rows += [json.loads(line) for line in fh if line.strip()]
    return rows


def fraction_label(f):
    return f if isinstance(f, str) else ("all" if f is None else f"{f:.0%}")


def oracle_summary(rows):
    groups = OrderedDict()
    for r in rows:
        groups.setdefault((r["group"], r["k"], fraction_label(r["fraction"])), []).append(r)
    out = []
    for (group, k, frac), rs in groups.items():
        solved = [r for r in rs if r["oracle"] is not None]
        out.append(OrderedDict(
            group=group, k=k, selection=frac, sets=len(rs),
            min_n=min(r["n"] for r in rs), max_n=max(r["n"] for r in rs),
            oracle_proven=len(solved),
            tamipami_matches_oracle=sum(r["match"] for r in solved),
            tamipami_worse_than_oracle=sum(tuple(r["tamipami"]) > tuple(r["oracle"]) for r in solved),
            tamipami_proven=sum(r["tamipami_proven"] for r in rs),
            candidates_limited=sum(r["candidates_limited"] for r in rs),
            max_tamipami_s=max(r["t_tamipami"] for r in rs),
        ))
    return out


def selection_summary(rows):
    by_key = OrderedDict()
    versions = list(OrderedDict.fromkeys(r["version"] for r in rows))
    for r in rows:
        by_key.setdefault((r["k"], fraction_label(r["fraction"])), {})[r["version"]] = r
    out = []
    for (k, frac), per_version in by_key.items():
        row = OrderedDict(k=k, selection=frac, n=next((r["n"] for r in per_version.values() if "n" in r), ""))
        for v in versions:
            r = per_version.get(v, {})
            if "error" in r:
                row[f"{v}: score"] = "failed: " + r["error"].split(":")[0]
                row[f"{v}: proven"] = row[f"{v}: s"] = row[f"{v}: MB"] = ""
                continue
            proven = r and not r["candidates_limited"] and all(s == "OPTIMAL" for s in r["stage_statuses"])
            row[f"{v}: score"] = f"{r['patterns']}/{r['degenerate_positions']}/{r['bases']}" if r else ""
            row[f"{v}: proven"] = proven if r else ""
            row[f"{v}: s"] = r.get("seconds", "")
            row[f"{v}: MB"] = r.get("peak_rss_mb", "")
        out.append(row)
    return out


def write(rows, path):
    if not rows:
        return
    fields = list(OrderedDict.fromkeys(f for r in rows for f in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\n{path}")
    widths = {f: max(len(f), *(len(str(r.get(f, ""))) for r in rows)) for f in fields}
    print("  ".join(f.ljust(widths[f]) for f in fields))
    for r in rows:
        print("  ".join(str(r.get(f, "")).ljust(widths[f]) for f in fields))


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "results"
    write(oracle_summary(load(os.path.join(d, "oracle_*.jsonl"))), os.path.join(d, "oracle_summary.csv"))
    write(selection_summary(load(os.path.join(d, "selection_*.jsonl"))), os.path.join(d, "selection_summary.csv"))
