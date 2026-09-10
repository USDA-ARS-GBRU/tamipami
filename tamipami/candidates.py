from __future__ import annotations
import logging
from typing import Dict, List, Set, Tuple, Optional
from itertools import product

import numpy as np
import scipy.cluster.hierarchy as sch

logger = logging.getLogger(__name__)

# --- IUPAC maps (complete for subsets of {A,C,G,T}) ---
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


# -------------------- Core helpers --------------------
def expand_iupac(code: str) -> Lis[str]:
    """Expand an IUPAC degenerate code into all concrete sequences."""
    pools: List[Tuple[str, ...]] = [IUPAC_TO_BASES[ch] for ch in code]
    seqs: List[str] = [""]
    for bases in pools:
        new_seqs: List[str] = []
        for p in seqs:
            for b in bases:
                new_seqs.append(p + b)
        seqs = new_seqs
    return seqs


def bases_to_iupac_symbol(bases: Set[str]) -> str:
    """Map a set of bases (size 1..4) to its IUPAC symbol."""
    key = tuple(sorted(bases))
    if key not in BASES_TO_IUPAC:
        raise ValueError(f"No IUPAC symbol for base set {bases}")
    return BASES_TO_IUPAC[key]


def unique_chars_per_position(strings: List[str]) -> List[Set[str]]:
    """Return per-position sets of bases from a list of equal-length sequences."""
    if not strings:
        return []
    k = len(strings[0])
    pos_sets: List[Set[str]] = [set() for _ in range(k)]
    for s in strings:
        for i, ch in enumerate(s):
            pos_sets[i].add(ch)
    return pos_sets


def expansion_is_subset_in_D(pos_sets: List[Set[str]], D_set: Set[str]) -> Optional[Set[str]]:
    """Return the expansion if it's subset of D_set; otherwise None (early exit)."""
    # Early reject if any position set is empty
    if any(len(ps) == 0 for ps in pos_sets):
        return None
    # Enumerate Cartesian product, short-circuit if a recombinant is missing
    exp: Set[str] = set()
    for tup in product(*[sorted(ps) for ps in pos_sets]):
        s = "".join(tup)
        if s not in D_set:
            return None
        exp.add(s)
    return exp


def code_for_pos_sets(pos_sets: List[Set[str]]) -> str:
    return "".join(bases_to_iupac_symbol(ps) for ps in pos_sets)


def degpos_complexity(code: str) -> Tuple[int, int]:
    """(#degenerate positions, total allowed bases) used for a 'best' code selector."""
    deg_positions = sum(1 for ch in code if len(IUPAC_TO_BASES[ch]) > 1)
    complexity = sum(len(IUPAC_TO_BASES[ch]) for ch in code)
    return deg_positions, complexity


# -------------------- Tree utilities --------------------
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
    if not strings:
        return None, None
    dm = condensed_hamming(strings)
    if dm.size == 0:
        return None, None
    Z = sch.linkage(dm, method=method)
    root = sch.to_tree(Z, rd=False)
    return Z, root


def collect_leaf_indices(node) -> List:
    """Collect leaf indices under a SciPy cluster node."""
    leaves: List[int] = []

    def _walk(n):
        if n.is_leaf():
            leaves.append(n.id)
        else:
            _walk(n.left)
            _walk(n.right)

    _walk(node)
    return leaves


def all_subtree_leaf_sets(root) -> List[List[int]]:
    """Enumerate leaves for every subtree (internal node + leaves)."""
    leaf_sets: List[List[int]] = []

    def _dfs(n):
        leaf_ids = collect_leaf_indices(n)
        leaf_sets.append(leaf_ids)
        if not n.is_leaf():
            _dfs(n.left)
            _dfs(n.right)

    _dfs(root)
    return leaf_sets


# -------------------- Candidate enumeration --------------------
def tree_induced_candidates(D: List[str], method: str = "average") -> List[Tuple[str, Set[str]]]:
    """
    Enumerate candidates from every subtree:
      - get leaf set L
      - compute positional sets S_i from L
      - code R = x_i S_i, accept iff expansion ⊆ D.
    """
    if not D:
        return []
    Z, root = build_dendrogram(D, method=method)
    if root is None:
        # Only singletons
        return [(s, {s}) for s in D]

    D_set = set(D)
    cands: List[Tuple[str, Set[str]]] = []
    for leaf_ids in all_subtree_leaf_sets(root):
        L = [D[i] for i in leaf_ids]
        pos_sets = unique_chars_per_position(L)
        exp = expansion_is_subset_in_D(pos_sets, D_set)
        if exp:
            code = code_for_pos_sets(pos_sets)
            cands.append((code, exp))
    return cands


def positional_subset_candidates(D: List[str]) -> List[Tuple[str, Set[str]]]:
    """
    Enumerate candidates by positional subsets:
      - For each position i, let cols[i] be bases observed in D
      - Enumerate all non-empty subsets of cols[i]
      - Combine across positions; accept iff expansion ⊆ D.
    """
    if not D:
        return []
    k = len(D[0])
    # Observed bases per position
    cols: List[Set[str]] = unique_chars_per_position(D)

    # Precompute non-empty subsets of each column's bases
    def nonempty_subsets(bs: Set[str]) -> List[Set[str]]:
        lst = sorted(bs)
        out: List[Set[str]] = []
        n = len(lst)
        for mask in range(1, 1 << n):
            out.append({lst[i] for i in range(n) if mask & (1 << i)})
        return out

    choices_per_pos: List[List[Set[str]]] = [nonempty_subsets(col) for col in cols]

    D_set = set(D)
    cands: List[Tuple[str, Set[str]]] = []

    # Enumerate product of positional choices
    for pos_choice in product(*choices_per_pos):
        exp = expansion_is_subset_in_D(list(pos_choice), D_set)
        if exp:
            code = code_for_pos_sets(list(pos_choice))
            cands.append((code, exp))
    return cands


def dedupe_candidates(cands: List[Tuple[str, Set[str]]]) -> List[Tuple[str, Set[str]]]:
    """
    Dedupe by expansion (frozenset). If multiple codes share the same expansion,
    keep the 'best' by (fewest degenerate positions, then fewest total bases).
    """
    by_exp: Dict[frozenset[str], Tuple[str, Set[str]]] = {}
    for code, exp in cands:
        key = frozenset(exp)
        best = by_exp.get(key)
        if best is None:
            by_exp[key] = (code, exp)
        else:
            # Choose better code
            (b_code, b_exp) = best
            b_score = degpos_complexity(b_code)
            c_score = degpos_complexity(code)
            if c_score < b_score:
                by_exp[key] = (code, exp)
    return list(by_exp.values())


def is_maximal_exact(code: str, D_set: Set[str]) -> bool:
    """
    Check if code is maximal exact wrt D_set:
      - for any position, adding any missing base should break exactness.
    """
    pos_sets = [set(IUPAC_TO_BASES[ch]) for ch in code]
    k = len(pos_sets)
    for i in range(k):
        missing = set(BASE_ORDER) - pos_sets[i]
        for b in missing:
            trial = [ps.copy() for ps in pos_sets]
            trial[i].add(b)
            if expansion_is_subset_in_D(trial, D_set):
                return False  # can grow while staying exact
    return True


def filter_maximal_exact(cands: List[Tuple[str, Set[str]]], D: List[str]) -> List[Tuple[str, Set[str]]]:
    D_set = set(D)
    out: List[Tuple[str, Set[str]]] = []
    for code, exp in cands:
        if is_maximal_exact(code, D_set):
            out.append((code, exp))
    return out


def enumerate_all_candidates(D: List[str], method: str = "average",
                             keep_only_maximal: bool = False) -> List[Tuple[str, Set[str]]]:
    """
    Full pipeline:
      1) Tree-induced rectangles
      2) Positional-subset rectangles
      3) Dedupe by expansion
      4) Optional: filter to maximal exact rectangles only
    Returns: list of (IUPAC code, expansion set) covering ALL exact rectangles.
    """
    tree_cands = tree_induced_candidates(D, method=method)
    pos_cands  = positional_subset_candidates(D)
    all_cands  = dedupe_candidates(tree_cands + pos_cands)
    if keep_only_maximal:
        all_cands = filter_maximal_exact(all_cands, D)
    return all_cands