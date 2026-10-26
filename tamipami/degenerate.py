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
import os
from typing import Dict, List, Set, Tuple, FrozenSet, Iterable, Optional

from ortools.sat.python import cp_model

logger = logging.getLogger(__name__)

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


def _default_num_workers() -> int:
    """CP-SAT worker count: TAMIPAMI_SOLVER_WORKERS if set, else up to 8 cores."""
    env = os.getenv("TAMIPAMI_SOLVER_WORKERS")
    if env:
        return max(1, int(env))
    return max(1, min(8, multiprocessing.cpu_count()))


def minimal_exact_cover(
    seqs: List[str],
    candidates_list: List[Tuple[str, Set[str]]],
    time_limit: float = 30.0,
    num_workers: Optional[int] = None,
    strict_exact: bool = False,  # Accept FEASIBLE when False
    deterministic_time: Optional[float] = None,
) -> List:
    """
    Exact-cover with lexicographic objective:

      Stage 1: Minimize number of selected patterns.
      Stage 2: Fix Stage 1 optimum; minimize total #degenerate positions.
      Stage 3: Fix Stage 1+2 optima; minimize total allowed bases.

    If strict_exact=False (default), return the best FEASIBLE assignment when OPTIMAL
    is not proven within 'time_limit'. Exact-cover constraints ensure every sequence
    is covered exactly once in all returned solutions.

    If `deterministic_time` is set, each stage is limited by CP-SAT's deterministic time
    (with interleaved search) instead of wall-clock time, so identical inputs give identical
    results on any machine; `time_limit` then only acts as a wall-clock backstop.
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
        deg_count = sum(1 for ch in code if len(IUPAC_CODES[ch]) > 1)
        deg_positions.append(deg_count)
        complexities.append(sum(len(IUPAC_CODES[ch]) for ch in code))

    # --- Build model & variables once; re-use across stages ---
    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x_{i}") for i in range(n_cands)]
    n_selected = cp_model.LinearExpr.Sum(x)
    total_degpos = cp_model.LinearExpr.WeightedSum(x, deg_positions)
    total_bases = cp_model.LinearExpr.WeightedSum(x, complexities)

    # Exact cover constraints: each sequence covered exactly once
    for seq_idx in range(n_seqs):
        model.AddExactlyOne(x[c] for c in covers[seq_idx])

    # Solver setup
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.relative_gap_limit = 0.0  # request optimality
    solver.parameters.num_search_workers = num_workers or _default_num_workers()
    solver.parameters.log_search_progress = False
    if deterministic_time is not None:
        solver.parameters.max_deterministic_time = deterministic_time
        solver.parameters.interleave_search = True

    def _extract_codes() -> List:
        return [candidates_list[i][0] for i in range(n_cands) if solver.Value(x[i]) == 1]

    # Lexicographic stages: (objective, constraint fixing its optimum for later stages)
    stages = [n_selected, total_degpos, total_bases]
    best_codes: Optional[List] = None
    for stage_no, objective in enumerate(stages, start=1):
        model.Minimize(objective)
        status = solver.Solve(model)
        logger.info("Stage %d status=%s objective=%s", stage_no, status, solver.ObjectiveValue())

        if status == cp_model.OPTIMAL:
            best_codes = _extract_codes()
            # Fix this optimum and start the next stage from the current solution,
            # which is feasible for it by construction
            model.Add(objective == int(round(solver.ObjectiveValue())))
            model.ClearHints()
            for var in x:
                model.AddHint(var, solver.Value(var))
        elif status == cp_model.FEASIBLE and not strict_exact:
            logger.warning(
                "Exact-cover accepted FEASIBLE at Stage %d (not proven OPTIMAL). "
                "Increase time_limit if strict optimality is required.", stage_no
            )
            return _extract_codes()
        elif best_codes is not None and not strict_exact:
            # A later stage found nothing in time; the previous stage's cover is still exact
            logger.warning(
                "Stage %d found no solution (status=%s); returning Stage %d result.",
                stage_no, status, stage_no - 1
            )
            return best_codes
        else:
            raise RuntimeError(
                f"Stage {stage_no} not optimal (status={status}). "
                "Increase time_limit or relax strict_exact."
            )
    return best_codes


# Deterministic CP-SAT budget per stage. On the RTW554 example, 20 was the smallest budget
# that found a solution for every selection tested; 30 adds margin (~6-7 s on an M-series Mac).
DETERMINISTIC_TIME = 30.0


def _cover_score(codes: List[str]) -> Tuple[int, int, int]:
    """Lexicographic objective of a cover: (#patterns, #degenerate positions, #bases)."""
    return (
        len(codes),
        sum(1 for code in codes for ch in code if len(IUPAC_CODES[ch]) > 1),
        sum(len(IUPAC_CODES[ch]) for code in codes for ch in code),
    )


# --------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------
def seqs_to_degenerates(
    seqs: Iterable[str],
    D: int = 2,                      # kept for API compatibility; ignored by enumeration
    time_limit: float = 60.0,
    deterministic_time: Optional[float] = DETERMINISTIC_TIME,
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
        time_limit: Wall-clock backstop per lexicographic stage (seconds). Results are only
            reproducible when `deterministic_time` is reached first.
        deterministic_time: CP-SAT deterministic time budget per stage; None uses
            wall-clock `time_limit` only (faster to tune, but not reproducible).
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
        # 1) Tree-induced candidates are few and always contain every singleton, so this
        #    solves quickly and gives a valid cover to fall back on and to hint the full solve.
        tree_cands = candidates.tree_induced_candidates(seqs, method=linkage)
        solve = dict(time_limit=time_limit, deterministic_time=deterministic_time)
        tree_codes = minimal_exact_cover(seqs=seqs, candidates_list=tree_cands, **solve)
        if len(tree_codes) == 1:
            return tree_codes  # a single exact rectangle is the unique optimum

        # 2) Full candidate set
        cand_list = candidates.combine_candidates(
            seqs, tree_cands, keep_only_maximal=keep_only_maximal
        )
        logger.info(
            "Enumerated %d candidates (maximal=%s, linkage=%s)",
            len(cand_list), keep_only_maximal, linkage
        )
        try:
            codes = minimal_exact_cover(seqs=seqs, candidates_list=cand_list, **solve)
        except RuntimeError as e:
            logger.warning("Full candidate solve failed (%s); using tree-induced cover.", e)
            return tree_codes
        # A time-limited FEASIBLE result is never allowed to be worse than the tree cover
        return min(codes, tree_codes, key=_cover_score)
    except Exception as e:
        logger.exception(f"An error occurred in seqs_to_degenerates: {e}")
        raise