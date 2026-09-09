#!/usr/bin/env python3
"""
Comprehensive tree-based benchmark (k=3..5) with old-style per-instance output.

Each instance:
  1) Samples reads of length k from 4^k (kept modest for speed).
  2) Builds ground-truth (GT) rectangles via dendrogram tree tiling:
     - Accept subtree rectangles only if expansion ⊆ D.
     - Guard against pairwise and sampled triplet merges.
     - Representatives (lexicographically smallest member) are Hamming≥2 apart.
     - Fill remainder with Hamming≥2 singletons if needed.
  3) Expands GT to dataset D and solves with tamipami.degenerate.seqs_to_degenerates.
  4) Evaluates coverage, ties, exact block matches, compactness sums, verdict.
  5) Writes per-instance JSON lines (sorted gt_codes & solver_codes) and summary CSV.

CLI:
  python benchmark_multi.py --k 3 --m 8 --samples 100 --outdir bench_multi_out
  python benchmark_multi.py --k 4 --m 5 --samples 100 --outdir bench_k4
  python benchmark_multi.py --k 5 --m 5 --samples 100 --outdir bench_k5 --linkage complete
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from itertools import product
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

import numpy as np
import scipy.cluster.hierarchy as sch

# Solver & candidates (your updated modules)
from tamipami.degenerate import seqs_to_degenerates
from tamipami import candidates as candmod

# ------------------ IUPAC maps ------------------
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
    # singletons
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
}
BASES_TO_IUPAC: Dict[Tuple[str, ...], str] = {tuple(sorted(v)): k for k, v in IUPAC_TO_BASES.items()}
BASE_ORDER: Tuple[str, ...] = ("A", "C", "G", "T")


# ------------------ helpers ------------------
def expand_iupac(code: str) -> List:
    """Expand an IUPAC degenerate code into all concrete sequences."""
    pools: List[Tuple[str, ...]] = [IUPAC_TO_BASES[ch] for ch in code]
    seqs: List[str] = [""]
    for bases in pools:
        new: List[str] = []
        for p in seqs:
            for b in bases:
                new.append(p + b)
        seqs = new
    return seqs


def bases_to_iupac_symbol(bases: Set[str]) -> str:
    """Map a set of bases (size 1..4) to its IUPAC symbol."""
    key = tuple(sorted(bases))
    if key not in BASES_TO_IUPAC:
        raise ValueError(f"No IUPAC symbol for base set {bases}")
    return BASES_TO_IUPAC[key]


def code_for_pos_sets(pos_sets: List[Set[str]]) -> str:
    """Build IUPAC code from per-position base sets."""
    return "".join(bases_to_iupac_symbol(ps) for ps in pos_sets)


def unique_chars_per_position(strings: List[str]) -> List[Set[str]]:
    """Return per-position sets of bases from equal-length sequences."""
    if not strings:
        return []
    k = len(strings[0])
    pos_sets: List[Set[str]] = [set() for _ in range(k)]
    for s in strings:
        for i, ch in enumerate(s):
            pos_sets[i].add(ch)
    return pos_sets


def expansion_is_subset_in_D(pos_sets: List[Set[str]], D_set: Set[str]) -> Optional[Set[str]]:
    """Return expansion if ⊆ D_set; else None (short-circuit on first missing recombinant)."""
    if any(len(ps) == 0 for ps in pos_sets):
        return None
    exp: Set[str] = set()
    for tup in product(*[sorted(ps) for ps in pos_sets]):
        s = "".join(tup)
        if s not in D_set:
            return None
        exp.add(s)
    return exp


def degpos_complexity(code: str) -> Tuple[int, int]:
    """Return (#degenerate positions, total allowed bases) for a code."""
    deg_positions = sum(1 for ch in code if len(IUPAC_TO_BASES[ch]) > 1)
    complexity = sum(len(IUPAC_TO_BASES[ch]) for ch in code)
    return deg_positions, complexity


def condensed_hamming(strings: List[str]) -> np.ndarray:
    """Build condensed Hamming distance matrix for SciPy linkage."""
    n = len(strings)
    dm = np.zeros(n * (n - 1) // 2, dtype=float)
    k = 0
    for i in range(n - 1):
        si = strings[i]
        for j in range(i + 1, n):
            sj = strings[j]
            dm[k] = sum(1 for a, b in zip(si, sj) if a != b)
            k += 1
    return dm


def build_dendrogram(strings: List[str], method: str = "average"):
    """Return linkage matrix and root node of the dendrogram."""
    dm = condensed_hamming(strings)
    Z = sch.linkage(dm, method=method)  # 'average' or 'complete' best for Hamming
    root = sch.to_tree(Z, rd=False)
    return Z, root


def collect_leaf_indices(node) -> List[int]:
    """Collect leaf indices under a SciPy cluster node."""
    leaves: List[int] = []

    def _walk(n) -> None:
        if n.is_leaf():
            leaves.append(n.id)
        else:
            _walk(n.left)
            _walk(n.right)

    _walk(node)
    return leaves


def is_disjoint_expansions(codes: List[str]) -> bool:
    """Check if expansions from codes are pairwise disjoint."""
    seen: Set[str] = set()
    for code in codes:
        for s in expand_iupac(code):
            if s in seen:
                return False
            seen.add(s)
    return True


def expansions_equal_setlist(gt_codes: List[str], other_codes: List[str]) -> bool:
    """
    Check if unordered sets of rectangle expansions are equal; enforce disjointness among other_codes expansions.
    """
    gt_sets = [frozenset(expand_iupac(c)) for c in gt_codes]
    other_sets = [frozenset(expand_iupac(c)) for c in other_codes]
    # enforce disjointness among other_sets
    for i in range(len(other_sets)):
        for j in range(i + 1, len(other_sets)):
            if other_sets[i].intersection(other_sets[j]):
                return False
    return set(other_sets) == set(gt_sets)


def union_is_exact(U: List[Set[str]], D_prime: Set[str]) -> bool:
    """True if Cartesian product of U expands entirely within D_prime."""
    return expansion_is_subset_in_D(U, D_prime) is not None


# ------------------ GT generator: tree tiling (generic k) ------------------
def exact_tree_tiling_gt(
    seqs: List[str],
    rng: random.Random,
    target_m: int,
    linkage: str = "average",
    min_rep_hamming: int = 2,
    check_triplets: bool = True,
) -> Tuple[List[str], List[int]]:
    """
    Build ground-truth rectangles from dendrogram subtrees:
      - accept subtree rectangle if expansion ⊆ D and anti-merge + disjointness guards hold
      - representatives (lexicographically smallest member) Hamming≥2 apart
      - fill remainder with Hamming≥2 singletons from D if needed
    Returns:
      gt_codes (IUPAC), exp_sizes (|expansion| for each GT code)
    """
    if not seqs:
        return [], []

    K = len(seqs[0])
    if any(len(s) != K for s in seqs):
        raise ValueError("All sequences must have the same length")

    D_set = set(seqs)
    _, root = build_dendrogram(seqs, method=linkage)

    queue = [root]
    gt_codes: List[str] = []
    exp_sizes: List[int] = []
    reps: List[str] = []         # lexicographically smallest sequence per accepted block
    covered_leaf_ids: Set[int] = set()
    used_set: Set[str] = set()   # tracks all sequences covered by accepted rectangles (disjointness)

    # -------------- BFS acceptance with guards --------------
    while queue and len(gt_codes) < target_m:
        node = queue.pop(0)
        leaf_ids = collect_leaf_indices(node)
        if all(i in covered_leaf_ids for i in leaf_ids):
            continue

        L = [seqs[i] for i in leaf_ids]
        pos_sets = unique_chars_per_position(L)
        cand_exp = expansion_is_subset_in_D(pos_sets, D_set)
        if cand_exp is None:
            # not exact -> split
            if not node.is_leaf():
                queue.append(node.left)
                queue.append(node.right)
            continue

        # --- Disjointness guard: reject if overlaps any previously accepted expansion ---
        if any(s in used_set for s in cand_exp):
            if not node.is_leaf():
                queue.append(node.left)
                queue.append(node.right)
            continue

        # --- Hamming-2 representatives guard ---
        rep = min(cand_exp)
        if any(sum(1 for a, b in zip(rep, r)) < min_rep_hamming for r in reps):
            if not node.is_leaf():
                queue.append(node.left)
                queue.append(node.right)
            continue

        # --- Pairwise anti-merge guard wrt already accepted rectangles ---
        D_prime = set().union(*[set(expand_iupac(c)) for c in gt_codes]) | cand_exp
        mergeable = False
        for code in gt_codes:
            ps_code = [set(IUPAC_TO_BASES[ch]) for ch in code]
            U = [ps_code[p] | pos_sets[p] for p in range(K)]
            if union_is_exact(U, D_prime):
                mergeable = True
                break
        if mergeable:
            if not node.is_leaf():
                queue.append(node.left)
                queue.append(node.right)
            continue

        # --- Optional sampled triplet anti-merge guard ---
        if check_triplets and len(gt_codes) >= 2:
            idxs = rng.sample(range(len(gt_codes)), k=min(6, len(gt_codes)))
            triplet_mergeable = False
            for i in range(0, len(idxs) - 1, 2):
                ci = gt_codes[idxs[i]]
                cj = gt_codes[idxs[i + 1]]
                ps_ci = [set(IUPAC_TO_BASES[ch]) for ch in ci]
                ps_cj = [set(IUPAC_TO_BASES[ch]) for ch in cj]
                U3 = [ps_ci[p] | ps_cj[p] | pos_sets[p] for p in range(K)]
                if union_is_exact(U3, D_prime):
                    triplet_mergeable = True
                    break
            if triplet_mergeable:
                if not node.is_leaf():
                    queue.append(node.left)
                    queue.append(node.right)
                continue

        # --- Accept rectangle ---
        code = code_for_pos_sets(pos_sets)
        gt_codes.append(code)
        exp_sizes.append(len(cand_exp))
        reps.append(rep)
        covered_leaf_ids |= set(leaf_ids)
        used_set |= cand_exp  # enforce disjointness

    # -------------- Hardened fallback: ensure len(gt_codes) == target_m --------------
    if len(gt_codes) < target_m:
        # Singletons available inside D that keep disjointness and rep H2
        remaining = target_m - len(gt_codes)
        # Prefer sequences in D that are not used yet
        candidates_singletons = [s for s in sorted(D_set) if s not in used_set]
        # Shuffle for variety
        rng.shuffle(candidates_singletons)

        for s in candidates_singletons:
            if remaining == 0:
                break
            # H2 vs existing reps
            if any(sum(1 for a, b in zip(s, r)) < min_rep_hamming for r in reps):
                continue
            # Disjointness trivially holds for singletons (sequence not in used_set)
            gt_codes.append(s)      # IUPAC singleton is the string itself
            exp_sizes.append(1)
            reps.append(s)
            used_set.add(s)
            remaining -= 1

        # As a last resort (should be rare), sample additional singletons from 4^K even if not in base reads,
        # but only if they are in D_set (we only cover D). This loop typically won't add anything new
        # if candidates_singletons exhausted; it's just belt-and-suspenders.
        if remaining > 0:
            universe = ["".join(p) for p in product(BASE_ORDER, repeat=K)]
            rng.shuffle(universe)
            for s in universe:
                if remaining == 0:
                    break
                if s not in D_set or s in used_set:
                    continue
                if any(sum(1 for a, b in zip(s, r)) < min_rep_hamming for r in reps):
                    continue
                gt_codes.append(s)
                exp_sizes.append(1)
                reps.append(s)
                used_set.add(s)
                remaining -= 1

        # Final assert: we must have reached target_m here; if not, raise clear error for diagnostics
        if len(gt_codes) < target_m:
            raise RuntimeError(
                f"GT tiler fallback failed: accepted={len(gt_codes)} < target_m={target_m}. "
                "Try reducing m or Hamming rep guard, or use signatures H2 GT mode."
            )

    return gt_codes, exp_sizes


# ------------------ per-instance runner ------------------
def run_instance_tree(
    example_idx: int,
    K: int,
    m: int,
    rng: random.Random,
    D_hamming: int,
    time_limit: float,
    linkage: str,
) -> Dict[str, object]:
    """
    One instance:
      - sample base reads of length K from 4^K,
      - build GT rectangles via tree tiling,
      - expand GT to D,
      - solve with seqs_to_degenerates,
      - evaluate and return old-style record (gt/solver sorted).
    """
    # Sample base reads (density scaled by K)
    universe = ["".join(p) for p in product(BASE_ORDER, repeat=K)]
    rng.shuffle(universe)
    size_ranges = {3: (24, 40), 4: (64, 120), 5: (96, 180)}
    low, high = size_ranges.get(K, (64, 160))
    base_reads = universe[:rng.randint(low, high)]

    # Build GT
    try:
        gt_codes, exp_sizes = exact_tree_tiling_gt(
            seqs=base_reads, rng=rng, target_m=m, linkage=linkage,
            min_rep_hamming=2, check_triplets=True,
        )
    except Exception as e:
        return {
            "example_idx": example_idx,
            "status": "failure",
            "mode": "rectangles",
            "k": K,
            "m": m,
            "gt_codes": [],
            "exp_sizes": [],
            "total_D_size": 0,
            "solver_codes": [],
            "solver_count": 0,
            "tie": False,
            "coverage_ok": False,
            "exact_block_matches": False,
            "verdict": "invalid",
            "gt_degpos_sum": 0,
            "gt_complexity_sum": 0,
            "solver_degpos_sum": 0,
            "solver_complexity_sum": 0,
            "error": f"GT generation error: {e!r}",
        }

    # Dataset D = union of expansions of GT codes
    D_set: Set[str] = set().union(*[set(expand_iupac(c)) for c in gt_codes])
    D: List[str] = sorted(D_set)
    total_D_size = len(D)

    # Solve with enumeration+CP-SAT
    try:
        solver_codes: List[str] = seqs_to_degenerates(
            D, D=D_hamming, time_limit=time_limit, linkage=linkage, keep_only_maximal=False
        )
    except Exception as e:
        return {
            "example_idx": example_idx,
            "status": "failure",
            "mode": "rectangles",
            "k": K,
            "m": m,
            "gt_codes": sorted(gt_codes),
            "exp_sizes": exp_sizes,
            "total_D_size": total_D_size,
            "solver_codes": [],
            "solver_count": 0,
            "tie": False,
            "coverage_ok": False,
            "exact_block_matches": False,
            "verdict": "invalid",
            "gt_degpos_sum": sum(degpos_complexity(c)[0] for c in gt_codes),
            "gt_complexity_sum": sum(degpos_complexity(c)[1] for c in gt_codes),
            "solver_degpos_sum": 0,
            "solver_complexity_sum": 0,
            "error": f"Solver exception: {e!r}",
        }

    # Evaluate coverage
    from tamipami.degenerate import minimal_exact_cover

    # 1) Enumerate candidates ONCE (same for solve and evaluation)
    all_cands = candmod.enumerate_all_candidates(
        D, method=linkage, keep_only_maximal=False
    )  # list[tuple[code:str, exp:set[str]]]
    exp_map: Dict[str, Set[str]] = {code: exp for code, exp in all_cands}

    # (Optional, Phase‑1 validation) Seed GT codes into the candidate pool
    # This guarantees the exact GT shapes are available to CP‑SAT among equal-expansion codes.
    existing_codes = set(exp_map.keys())
    for code in gt_codes:
        if code not in existing_codes:
            exp = set(expand_iupac(code))
            if exp.issubset(set(D)):
                all_cands.append((code, exp))
                exp_map[code] = exp
                existing_codes.add(code)

    # 2) Solve lexicographically, OPTIMAL-only, using the SAME candidate list
    try:
        solver_codes: List[str] = minimal_exact_cover(
            seqs=D, candidates_list=all_cands, time_limit=time_limit
        )
    except Exception as e:
        return {
            "example_idx": example_idx, "status": "failure",
            # ... (same failure record fields you already return) ...
            "error": f"Solver exception: {e!r}",
        }

    # 3) Evaluate coverage using exp_map built from all_cands (same universe)
    solver_sets = [exp_map.get(code, set(expand_iupac(code))) for code in solver_codes]
    disjoint_solver = True
    for i in range(len(solver_sets)):
        for j in range(i + 1, len(solver_sets)):
            if solver_sets[i].intersection(solver_sets[j]):
                disjoint_solver = False
                break
        if not disjoint_solver:
            break
    union_solver: Set[str] = set().union(*solver_sets) if solver_sets else set()
    coverage_ok = (disjoint_solver and union_solver == set(D))

    # 4) Diagnostics: candidate_count, gt_missing_in_candidates (by code presence)
    candidate_count = len(all_cands)
    gt_missing_in_candidates = sorted([c for c in gt_codes if c not in existing_codes])

    # Sort codes for easier comparison
    gt_codes_sorted = sorted(gt_codes)
    solver_codes_sorted = sorted(solver_codes)

    solver_count = len(solver_codes_sorted)
    tie = (solver_count == len(gt_codes_sorted))
    exact_block_matches = expansions_equal_setlist(gt_codes_sorted, solver_codes_sorted) if coverage_ok and tie else False

    # --- Diagnostics (Patch A) ---
    # Are any GT codes missing from the candidate pool?
    cand_code_set = set(exp_map.keys())
    missing_gt = [c for c in gt_codes_sorted if c not in cand_code_set]

    # Optional k=3 oracle: true minimal cardinality using complete IUPAC enumeration
    oracle_min = -1
    oracle_status = "disabled"
    try:
        if K == 3:
            # You must provide oracle_min_cardinality_k3(D) in the file or import it.
            from tamipami.oracle_k3 import oracle_min_cardinality_k3  # <— or local function
            oracle_min, oracle_status = oracle_min_cardinality_k3(D, time_limit=10.0)
    except Exception as _e:
        oracle_min, oracle_status = -1, f"error:{_e!r}"

    # Verdict
    if coverage_ok:
        if tie and exact_block_matches:
            verdict = "equivalent"
        elif solver_count < len(gt_codes_sorted):
            verdict = "better"
        elif solver_count > len(gt_codes_sorted):
            # If oracle exists and proves fewer rectangles are possible, flag candidate-gap
            if oracle_min != -1 and solver_count > oracle_min:
                verdict = "worse_candidate_gap"
            else:
                verdict = "worse"
        else:
            verdict = "valid-different"
    else:
        verdict = "invalid"

    # Compactness sums
    gt_degpos_sum = sum(degpos_complexity(c)[0] for c in gt_codes_sorted)
    gt_complexity_sum = sum(degpos_complexity(c)[1] for c in gt_codes_sorted)
    solver_degpos_sum = sum(degpos_complexity(c)[0] for c in solver_codes_sorted)
    solver_complexity_sum = sum(degpos_complexity(c)[1] for c in solver_codes_sorted)

    return {
        "example_idx": example_idx,
        "status": "success",
        "mode": "rectangles",
        "k": K,
        "m": m,
        "gt_codes": gt_codes_sorted,             # sorted
        "exp_sizes": exp_sizes,
        "total_D_size": total_D_size,
        "solver_codes": solver_codes_sorted,     # sorted
        "solver_count": solver_count,
        "tie": tie,
        "coverage_ok": coverage_ok,
        "exact_block_matches": exact_block_matches,
        "verdict": verdict,
        "gt_degpos_sum": gt_degpos_sum,
        "gt_complexity_sum": gt_complexity_sum,
        "solver_degpos_sum": solver_degpos_sum,
        "solver_complexity_sum": solver_complexity_sum,
        # --- Diagnostics in the returned record ---
        "gt_missing_in_candidates": sorted(missing_gt),
        "candidate_count": len(all_cands),
        "oracle_min": oracle_min,
        "oracle_status": oracle_status,
        "error": "",
    }


# ------------------ I/O ------------------
def write_jsonl(path: Path, records: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")


def write_summary_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ------------------ CLI ------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Comprehensive tamipami benchmark (old-style output, k=3..5)")
    parser.add_argument("--k", type=int, default=3, help="Sequence length k (3..5).")
    parser.add_argument("--m", type=int, default=8, help="Target GT rectangle count per instance.")
    parser.add_argument("--samples", type=int, default=100, help="Number of instances to run.")
    parser.add_argument("--D", type=int, default=1, help="Hamming radius (kept for API compatibility).")
    parser.add_argument("--time-limit", type=float, default=30.0, help="CP-SAT time limit per instance.")
    parser.add_argument("--linkage", type=str, default="average", choices=["average", "complete"],
                        help="Dendrogram linkage used during candidate enumeration (tree-induced path).")
    parser.add_argument("--outdir", type=Path, default=Path("bench_multi_out"),
                        help="Output directory.")
    args = parser.parse_args()

    rng_master = random.Random(12345)
    print(f"\n=== Tree-based benchmark for k={args.k}, m={args.m}, samples={args.samples}, ===")

    per_instances: List[Dict[str, object]] = []
    successes = failures = 0
    ties = coverage_ok = exact_block_matches = 0
    verdict_counts = {
        "equivalent": 0, "better": 0, "worse": 0, "valid-different": 0, "invalid": 0,
    }

    for i in range(args.samples):
        rng = random.Random(rng_master.getrandbits(64))
        rec = run_instance_tree(
            example_idx=i,
            K=args.k,
            m=args.m,
            rng=rng,
            D_hamming=args.D,
            time_limit=args.time_limit,
            linkage=args.linkage,
        )
        per_instances.append(rec)

        if rec["status"] == "success":
            successes += 1
            ties += int(rec["tie"])                 # type: ignore[arg-type]
            coverage_ok += int(rec["coverage_ok"])  # type: ignore[arg-type]
            exact_block_matches += int(rec["exact_block_matches"])  # type: ignore[arg-type]
            v = str(rec.get("verdict", "invalid"))
            if v in verdict_counts:
                verdict_counts[v] += 1
        else:
            failures += 1

        if (i + 1) % max(1, args.samples // 10) == 0:
            print(f"  processed {i+1}/{args.samples} ...")

    outdir_k = args.outdir / f"k{args.k}_m{args.m}"
    write_jsonl(outdir_k / "instances.jsonl", per_instances)

    # Old-style summary row
    summary_row: Dict[str, object] = {
        "k": args.k,
        "m": args.m,
        "samples": args.samples,
        "successes": successes,
        "failures": failures,
        "ties": ties,
        "coverage_ok": coverage_ok,
        "exact_block_matches": exact_block_matches,
        "equivalent": verdict_counts["equivalent"],
        "better": verdict_counts["better"],
        "worse": verdict_counts["worse"],
        "valid-different": verdict_counts["valid-different"],
        "invalid": verdict_counts["invalid"],
    }
    write_summary_csv(args.outdir / "summary.csv", [summary_row])

    print(
        f"Summary k={args.k}: successes={successes}/{args.samples}, ties={ties}, "
        f"coverage_ok={coverage_ok}, exact_block_matches={exact_block_matches}, failures={failures}"
    )
    print(
        f"Breakdown ⇒ equivalent={verdict_counts['equivalent']}, better={verdict_counts['better']}, "
        f"worse={verdict_counts['worse']}, valid-different={verdict_counts['valid-different']}, invalid={verdict_counts['invalid']}"
    )
    print("\nWrote outputs to:", args.outdir.resolve())


if __name__ == "__main__":
    main()