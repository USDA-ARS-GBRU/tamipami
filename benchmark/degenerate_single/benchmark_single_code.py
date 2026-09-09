#!/usr/bin/env python3
"""
Empirical benchmarking: single-degenerate-code case.

For k in {3,4,5}, sample N random degenerate IUPAC codes (each with at least one
degenerate position), expand to the full set of sequences, run the solver
(tamipami.degenerate.seqs_to_degenerates), and record successes/failures.

Success criteria (single-degenerate benchmark):
  - Solver returns exactly ONE code, and
  - Expanding that solver code yields exactly the same set as the input expansion.

Outputs (per k):
  - sampled_iupac_k{k}.txt: full list of sampled degenerate input codes
  - expanded_sets_k{k}.jsonl: JSON lines; each line has {"code": "...", "seqs": [...]}
  - solver_outputs_unique_k{k}.txt: unique solver codes seen
  - failures_k{k}.csv: failure rows with:
        input_code, k, expanded_size,
        solver_output_codes,
        solver_union_size,
        mismatch_details,
        expanded_set (input expansion),
        solver_union_expanded_set (union of all solver output expansions)
  - summary.json: successes/failures per k

Run:
  python benchmark_single_code.py --samples 1000 --seed 42 --D 2 --time-limit 30 --outdir bench_out
"""

from __future__ import annotations
import argparse
import csv
import json
import random
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

# --- Import the solver from tamipami.degenerate ---
_import_errors: List[str] = []
try:
    from tamipami.degenerate import seqs_to_degenerates  # type: ignore
except Exception as e:
    _import_errors.append(f"tamipami.degenerate import failed: {e}")
    seqs_to_degenerates = None  # type: ignore


# --- IUPAC map (complete for sizes 1..4) ---
IUPAC_TO_BASES: Dict[str, Tuple[str, ...]] = {
    # 3-letter sets
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    # 2-letter sets
    "K": ("G", "T"),
    "M": ("A", "C"),
    "R": ("A", "G"),
    "S": ("C", "G"),
    "W": ("A", "T"),
    "Y": ("C", "T"),
    # 4-letter set
    "N": ("A", "C", "G", "T"),
    # 1-letter sets
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
}

DEGENERATE_SYMBOLS: Tuple[str, ...] = tuple(
    s for s, bases in IUPAC_TO_BASES.items() if len(bases) > 1
)
NON_DEGENERATE_SYMBOLS: Tuple[str, ...] = ("A", "C", "G", "T")

# Weights to get a reasonable spread of degeneracy types
DEGENERATE_WEIGHTS: Dict[str, float] = {
    # 3-letter
    "B": 1.0, "D": 1.0, "H": 1.0, "V": 1.0,
    # 2-letter
    "K": 2.0, "M": 2.0, "R": 2.0, "S": 2.0, "W": 2.0, "Y": 2.0,
    # 4-letter
    "N": 1.5,
}
NON_DEGENERATE_WEIGHT: float = 1.0  # equal for A/C/G/T


def expand_iupac(code: str) -> List:
    """Expand an IUPAC degenerate code into all concrete sequences."""
    pools: List[Tuple[str, ...]] = [IUPAC_TO_BASES[ch] for ch in code]
    seqs = [""]
    for pos_bases in pools:
        new_seqs: List[str] = []
        for prefix in seqs:
            for b in pos_bases:
                new_seqs.append(prefix + b)
        seqs = new_seqs
    return seqs


def is_degenerate(code: str) -> bool:
    """True if at least one position is degenerate (allowed set size > 1)."""
    return any(len(IUPAC_TO_BASES[ch]) > 1 for ch in code)


def sample_random_iupac(k: int, rng: random.Random, p_degenerate: float = 0.6) -> str:
    """
    Sample an IUPAC code of length k.
    Ensures at least one degenerate position.

    p_degenerate: per-position probability to choose from degenerate symbols.
    """
    code_chars: List[str] = []
    for _ in range(k):
        if rng.random() < p_degenerate:
            symbols = list(DEGENERATE_WEIGHTS.keys())
            weights = [DEGENERATE_WEIGHTS[s] for s in symbols]
            ch = rng.choices(symbols, weights=weights, k=1)[0]
        else:
            ch = rng.choice(NON_DEGENERATE_SYMBOLS)
        code_chars.append(ch)

    # Ensure at least one degenerate position
    if not any(len(IUPAC_TO_BASES[ch]) > 1 for ch in code_chars):
        idx = rng.randrange(k)
        symbols = list(DEGENERATE_WEIGHTS.keys())
        weights = [DEGENERATE_WEIGHTS[s] for s in symbols]
        code_chars[idx] = rng.choices(symbols, weights=weights, k=1)[0]

    return "".join(code_chars)


def solver_available() -> bool:
    return callable(seqs_to_degenerates)


def run_case(code: str, D_hamming: int, time_limit: float) -> Dict:
    """
    Run a single benchmark case:
      - Expand code to sequences S_in
      - Run solver on S_in
      - Check success criteria (1 code whose expansion exactly equals S_in)

    Returns a record dict with fields:
      input_code, k, expanded_size,
      solver_output_codes, success, mismatch_details,
      solver_union_size
    """
    S_in = expand_iupac(code)
    k = len(code)

    # Execute solver
    try:
        solver_codes: List[str] = seqs_to_degenerates(S_in, D=D_hamming, time_limit=time_limit)
    except Exception as e:
        return {
            "input_code": code,
            "k": k,
            "expanded_size": len(S_in),
            "solver_output_codes": [],
            "success": False,
            "mismatch_details": f"Solver raised exception: {e!r}",
            "solver_union_size": 0,
        }

    # Success criteria:
    # (1) Exactly one solver code
    # (2) Expanded solver code equals set(S_in)
    if len(solver_codes) == 1:
        out_code = solver_codes[0]
        S_out = set(expand_iupac(out_code))
        S_in_set = set(S_in)
        if S_out == S_in_set:
            return {
                "input_code": code,
                "k": k,
                "expanded_size": len(S_in),
                "solver_output_codes": solver_codes,
                "success": True,
                "mismatch_details": "",
                "solver_union_size": len(S_out),
            }
        else:
            return {
                "input_code": code,
                "k": k,
                "expanded_size": len(S_in),
                "solver_output_codes": solver_codes,
                "success": False,
                "mismatch_details": (
                    "Expanded solver code differs from input expansion "
                    f"(input_size={len(S_in_set)}, output_size={len(S_out)})"
                ),
                "solver_union_size": len(S_out),
            }
    else:
        # Build union of solver expansions (for inspection)
        union: Set[str] = set()
        for c in solver_codes:
            union.update(expand_iupac(c))
        return {
            "input_code": code,
            "k": k,
            "expanded_size": len(S_in),
            "solver_output_codes": solver_codes,
            "success": False,
            "mismatch_details": f"Expected 1 code; got {len(solver_codes)}",
            "solver_union_size": len(union),
        }


def write_list(path: Path, items: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for item in items:
            f.write(f"{item}\n")


def write_expanded_jsonl(path: Path, records: Sequence[Tuple[str, List[str]]]) -> None:
    """
    Write JSONL: one object per line with {"code": "<IUPAC>", "seqs": [ ... ]}.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for code, seqs in records:
            obj = {"code": code, "seqs": seqs}
            f.write(json.dumps(obj) + "\n")


def write_failures_csv(path: Path, failures: List[Dict]) -> None:
    """
    Failure CSV fields:
      input_code, k, expanded_size,
      solver_output_codes, solver_union_size,
      mismatch_details,
      expanded_set, solver_union_expanded_set
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "input_code",
        "k",
        "expanded_size",
        "solver_output_codes",
        "solver_union_size",
        "mismatch_details",
        "expanded_set",
        "solver_union_expanded_set",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for rec in failures:
            input_code = rec["input_code"]
            expanded_set = expand_iupac(input_code)
            solver_output_codes = rec.get("solver_output_codes", [])
            # Union of solver expansions (if any)
            union_exp: Set[str] = set()
            for c in solver_output_codes:
                union_exp.update(expand_iupac(c))

            w.writerow({
                "input_code": input_code,
                "k": rec["k"],
                "expanded_size": rec["expanded_size"],
                "solver_output_codes": "|".join(solver_output_codes),
                "solver_union_size": rec.get("solver_union_size", 0),
                "mismatch_details": rec.get("mismatch_details", ""),
                "expanded_set": " ".join(expanded_set),
                "solver_union_expanded_set": " ".join(sorted(union_exp)),
            })


def main():
    parser = argparse.ArgumentParser(description="Benchmark single-degenerate-code cases.")
    parser.add_argument("--samples", type=int, default=1000, help="Samples per k (3,4,5).")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed.")
    parser.add_argument("--D", type=int, default=2, help="Hamming radius used by candidate generation (D).")
    parser.add_argument("--time-limit", type=float, default=30.0, help="CP-SAT time limit per case (seconds).")
    parser.add_argument("--outdir", type=Path, default=Path("benchmark_single_code_out"), help="Output directory.")
    args = parser.parse_args()

    rng = random.Random(args.seed)

    if not solver_available():
        raise RuntimeError(
            "Could not import tamipami.degenerate.seqs_to_degenerates. "
            f"Import errors: {', '.join(_import_errors)}"
        )

    summary: Dict[int, Dict[str, int]] = {}

    for k in (3, 4, 5):
        print(f"\n=== Benchmarking k={k} with {args.samples} random degenerate codes ===")

        sampled_codes: List[str] = []
        expanded_records: List[Tuple[str, List[str]]] = []

        successes = 0
        failures: List[Dict] = []
        unique_solver_outputs: Set[str] = set()

        for i in range(args.samples):
            code = sample_random_iupac(k, rng)
            sampled_codes.append(code)

            # Store expanded set for "full set of degenerate sequences generated"
            expanded_seqs = expand_iupac(code)
            expanded_records.append((code, expanded_seqs))

            rec = run_case(code, D_hamming=args.D, time_limit=args.time_limit)

            if rec["success"]:
                successes += 1
            else:
                failures.append(rec)

            for c in rec["solver_output_codes"]:
                unique_solver_outputs.add(c)

            if (i + 1) % 100 == 0:
                print(f"  processed {i+1}/{args.samples} ...")

        # Write sampled input codes (full list per k)
        write_list(args.outdir / f"sampled_iupac_k{k}.txt", sampled_codes)

        # Write expanded sets for EVERY sampled code (JSONL)
        write_expanded_jsonl(args.outdir / f"expanded_sets_k{k}.jsonl", expanded_records)

        # Write unique solver outputs per k
        write_list(args.outdir / f"solver_outputs_unique_k{k}.txt", sorted(unique_solver_outputs))

        # Write failures CSV with details and expanded sets + solver union
        if failures:
            write_failures_csv(args.outdir / f"failures_k{k}.csv", failures)

        # Per-k summary
        summary[k] = {
            "samples": args.samples,
            "successes": successes,
            "failures": len(failures),
        }
        print(f"Summary k={k}: successes={successes}, failures={len(failures)}")

    # Write overall summary JSON
    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print("\nWrote outputs to:", args.outdir.resolve())


if __name__ == "__main__":
    main()