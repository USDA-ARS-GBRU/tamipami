"""
tamipami.degenerate
-------------------

Exact-cover solver for degenerate IUPAC rectangles.

- Candidate generation is delegated to tamipami.candidates.enumerate_all_candidates()
  (tree-induced + positional-subsets + dedupe; optional maximal filter).
- CP-SAT solve uses a three-stage lexicographic objective and requires OPTIMAL solutions:
    1) Minimize number of rectangles
    2) Minimize total degenerate positions
    3) Minimize total allowed bases

Notes:
- The legacy rectangle-expansion / pair-seeding code has been removed.
- The 'D' parameter in seqs_to_degenerates is retained for API compatibility
  but ignored, since exact enumeration does not rely on Hamming radius.
"""

from __future__ import annotations

import logging
import multiprocessing
from typing import Dict, List, Set, Tuple, FrozenSet, Iterable, Optional

from ortools.sat.python import cp_model

# --- Import exact candidate enumerator ---
try:
    # Package-relative import (preferred)
    from tamipami import candidates
except ImportError:  # pragma: no cover
    # Fallback if executed outside package context
    import candidates

# --- Minimal IUPAC map for solver metrics (degenerate position counts, complexity) ---
IUPAC_CODES: Dict[str, Set[str]] = {
    # 3-letter sets
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    # 2-letter sets
    "K": {"G", "T"},
    "M": {"A", "C"},
    "R": {"A", "G"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "Y": {"C", "T"},
    # 4-letter set
    "N": {"A", "C", "G", "T"},
    # singletons
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
}
IUPAC_LOOKUP: Dict[FrozenSet[str], str] = {frozenset(v): k for k, v in IUPAC_CODES.items()}


# --------------------------------------------------------------------------------------
# Lexicographic, OPTIMAL-only exact-cover solver
# --------------------------------------------------------------------------------------
def minimal_exact_cover(
    seqs: List[str],
    candidates_list: List[Tuple[str, Set[str]]],
    time_limit: float = 30.0,
    num_workers: Optional[int] = None,
) -> List:
    """
    Exact-cover with lexicographic objective (OPTIMAL-only):

    Stage 1: Minimize number of selected patterns
    Stage 2: Fix Stage 1 optimum; minimize total #degenerate positions across selected patterns
    Stage 3: Fix Stage 1+2 optima; minimize total allowed bases across all positions

    Raises if any stage is not OPTIMAL (feasible-but-not-proven-minimal solutions are rejected).

    Args:
        seqs: Sorted list of sequences to cover (order doesn't matter; uniqueness does).
        candidates_list: List of (code, covered_set) where:
            - code: IUPAC degenerate string (e.g., "MRT...")
            - covered_set: set of sequences this code expands to and covers
        time_limit: CP-SAT solver time per stage (seconds).
        num_workers: optional worker threads; defaults to CPU count.

    Returns:
        Listchosen candidate codes forming the lexicographically minimal exact cover.
    """
    # --- Indexing & incidence ---
    seq_index = {s: i for i, s in enumerate(seqs)}
    n_seqs = len(seqs)
    n_cands = len(candidates_list)

    if n_seqs == 0:
        return []

    covers: List[List[int]] = [[] for _ in range(n_seqs)]
    for cand_idx, (_, covered_set) in enumerate(candidates_list):
        for s in covered_set:
            idx = seq_index.get(s)
            if idx is not None:
                covers[idx].append(cand_idx)

    # Defensive: every sequence must be coverable by at least one candidate
    uncovered = [seqs[i] for i, lst in enumerate(covers) if not lst]
    if uncovered:
        raise ValueError(
            f"Infeasible: {len(uncovered)} sequences are uncovered by any candidate; "
            f"example: {uncovered[:3]}"
        )

    # --- Per-candidate metrics ---
    deg_positions: List[int] = []
    complexities: List[int] = []

    for code, _ in candidates_list:
        # number of positions that are degenerate (allowed set size > 1)
        deg_count = sum(1 for ch in code if len(IUPAC_CODES[ch]) > 1)
        deg_positions.append(deg_count)
        
        # total allowed bases across all positions (sum of set sizes)
        complexities.append(sum(len(IUPAC_CODES[ch]) for ch in code))


    # --- Build model & variables once; re-use across stages ---
    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x_{i}") for i in range(n_cands)]

    # Exact cover constraints: each sequence covered exactly once
    for seq_idx in range(n_seqs):
        model.Add(sum(x[c] for c in covers[seq_idx]) == 1)

    # Solver setup
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.relative_gap_limit = 0.0  # require exact optimality
    solver.parameters.num_search_workers = num_workers or multiprocessing.cpu_count()
    solver.parameters.log_search_progress = False

    # --- Stage 1: minimize number of selected patterns ---
    model.Minimize(sum(x))
    status = solver.Solve(model)
    logging.info("Stage 1 status=%s objective=%s", status, solver.ObjectiveValue())

    if status != cp_model.OPTIMAL:
        raise RuntimeError(
            f"Stage 1 not optimal (status={status}). Increase time_limit or reduce candidates."
        )
    min_patterns = int(solver.ObjectiveValue())

    # --- Stage 2: fix #patterns; minimize total degenerate positions ---
    model.Add(sum(x) == min_patterns)
    model.Minimize(sum(x[i] * deg_positions[i] for i in range(n_cands)))
    status = solver.Solve(model)
    logging.info("Stage 2 status=%s objective=%s", status, solver.ObjectiveValue())

    if status != cp_model.OPTIMAL:
        raise RuntimeError(
            f"Stage 2 not optimal (status={status}). Increase time_limit or reduce candidates."
        )
    min_degpos = int(solver.ObjectiveValue())

    # --- Stage 3: fix previous minima; minimize total bases ---
    model.Add(sum(x[i] * deg_positions[i] for i in range(n_cands)) == min_degpos)
    model.Minimize(sum(x[i] * complexities[i] for i in range(n_cands)))
    status = solver.Solve(model)
    logging.info("Stage 3 status=%s objective=%s", status, solver.ObjectiveValue())

    if status != cp_model.OPTIMAL:
        raise RuntimeError(
            f"Stage 3 not optimal (status={status}). Increase time_limit or reduce candidates."
        )

    # Extract chosen patterns
    chosen_codes = [candidates_list[i][0] for i in range(n_cands) if solver.Value(x[i]) == 1]
    return chosen_codes


# --------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------
def seqs_to_degenerates(
    seqs: Iterable[str],
    D: int = 2,                      # kept for API compatibility; ignored by enumeration
    time_limit: float = 30.0,
    linkage: str = "average",        # 'average' (default) or 'complete' preferred for Hamming
    keep_only_maximal: bool = False  # speed heuristic; default False to preserve optimality options
) -> List:
    """
    Build a (near-)complete candidate set from the dataset and solve exact cover with CP-SAT.

    Candidate generation uses tamipami.candidates.enumerate_all_candidates(...)
    (tree-induced + positional-subsets + dedupe; optional maximal filter).

    Args:
        seqs: Equal-length sequences over A,C,G,T.
        D: Ignored (prefilter replaced by exact enumeration); kept for API compatibility.
        time_limit: CP-SAT time per lexicographic stage (seconds).
        linkage: Dendrogram linkage for tree-induced candidates ('average' or 'complete').
        keep_only_maximal: If True, restrict candidates to maximal exact rectangles (heuristic).

    Returns:
        ListIUPAC codes forming an exact cover of the dataset, lexicographically optimal.
    """
    seqs = sorted(set(seqs))
    if not seqs:
        return []
    if len({len(s) for s in seqs}) != 1:
        raise ValueError("All sequences must have the same length")
    if len(seqs) == 1:
        return seqs

    try:
        cand_list = candidates.enumerate_all_candidates(
            seqs, method=linkage, keep_only_maximal=keep_only_maximal
        )
        logging.info(
            "Enumerated %d candidates (maximal=%s, linkage=%s)",
            len(cand_list), keep_only_maximal, linkage
        )
        return minimal_exact_cover(seqs=seqs, candidates_list=cand_list, time_limit=time_limit)
    except Exception as e:
        logging.exception(f"An error occurred in seqs_to_degenerates: {e}")
        raise