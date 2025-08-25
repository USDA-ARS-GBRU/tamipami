import logging
from typing import Dict, List, Set, Tuple, FrozenSet, Iterable, Optional

from ortools.sat.python import cp_model

# Package version (resolved dynamically when installed; fallback for local use)
try:
    from importlib.metadata import version, PackageNotFoundError  # type: ignore
    try:
        __version__ = version("exact-degenerate")
    except PackageNotFoundError:
        __version__ = "0.0.0"
except Exception:  # pragma: no cover - very defensive
    __version__ = "0.0.0"

# IUPAC code map (includes 1-, 2-, 3-, and 4-letter sets)
IUPAC_CODES: Dict[str, Set[str]] = {
    "B": {"G", "T", "C"},
    "D": {"G", "A", "T"},
    "H": {"A", "T", "C"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "N": {"G", "A", "T", "C"},
    "R": {"G", "A"},
    "S": {"G", "C"},
    "V": {"G", "A", "C"},
    "W": {"A", "T"},
    "Y": {"T", "C"},
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
}
IUPAC_LOOKUP: Dict[FrozenSet[str], str] = {
    frozenset(v): k for k, v in IUPAC_CODES.items()
}

BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}
INT_TO_BASE = {v: k for k, v in BASE_TO_INT.items()}


def encode_seq(s: str, pow4: List[int]) -> int:
    """Encode a DNA sequence as a base-4 integer.

    Maps A,C,G,T to 0,1,2,3 and sums BASE_TO_INT[ch] * pow4[i] for each position i
    (starting at 0). Typically, pow4[i] equals 4**i.

    Args:
        s (str): DNA sequence to encode.
        pow4 (List[int]): Precomputed powers of 4; length must be at least len(s).

    Returns:
        int: Base-4 integer index of the sequence.

    Raises:
        KeyError: If s contains characters other than A,C,G,T.
        IndexError: If pow4 is shorter than len(s).
    """
    x = 0
    for i, ch in enumerate(s):
        x += BASE_TO_INT[ch] * pow4[i]
    return x


def decode_seq(x: int, L: int, pow4: List[int]) -> str:
    """
    Decode a base-4 encoded integer into a DNA sequence.

    Args:
        x: Non-negative integer representing the sequence in base-4.
        L: Length of the decoded sequence.
        pow4: Precomputed powers of 4 where pow4[i] == 4**i and len(pow4) >= L.

    Returns:
        str: DNA string of length L decoded using the INT_TO_BASE mapping.
    """
    chars = []
    for i in range(L):
        d = (x // pow4[i]) % 4
        chars.append(INT_TO_BASE[d])
    return "".join(chars)


def neighbors_within_hamming(
    idx: int, L: int, pow4: List[int], digits: List[int], D: int
) -> List[int]:
    """
    Enumerate neighbor indices within Hamming distance up to D for a base-4, 
    length-L encoding using single- and, if requested, double-position digit flips.

    Args:
        idx: Base-4 encoded index of the current sequence.
        L: Sequence length (number of base-4 digits).
        pow4: Precomputed powers of 4 where pow4[p] == 4**p.
        digits: Base-4 digits of idx for positions 0..L-1 (each in {0, 1, 2, 3}).
        D: Maximum Hamming distance to enumerate; only 0, 1, or 2 are considered (values >2 are treated as 2).

    Returns:
        List[int]: Neighbor indices at distance 1 (and 2 if D >= 2). The original index is excluded; 
        returns an empty list when D == 0. Order is deterministic but not guaranteed to be sorted.
    """
    nbs: List[int] = []
    # Precompute single-position alternatives
    single = []
    for pos in range(L):
        cur = digits[pos]
        base = idx - cur * pow4[pos]
        for alt in (0, 1, 2, 3):
            if alt == cur:
                continue
            single.append(base + alt * pow4[pos])
    if D >= 1:
        nbs.extend(single)
    if D >= 2:
        # Double flips (pairs of positions)
        for a in range(L):
            cur_a = digits[a]
            base_a = idx - cur_a * pow4[a]
            for alt_a in (0, 1, 2, 3):
                if alt_a == cur_a:
                    continue
                idx_a = base_a + alt_a * pow4[a]
                for b in range(a + 1, L):
                    cur_b = digits[b]
                    base_b = idx_a - cur_b * pow4[b]
                    for alt_b in (0, 1, 2, 3):
                        if alt_b == cur_b:
                            continue
                        nbs.append(base_b + alt_b * pow4[b])
    return nbs


def iupac_from_pos_sets(pos_sets: Tuple[FrozenSet[str], ...]) -> Optional[str]:
    """
    Convert per-position base sets into an IUPAC-degenerate sequence.

    Each frozen set in pos_sets is mapped to its IUPAC symbol via IUPAC_LOOKUP
    and the symbols are concatenated. If any set lacks a corresponding symbol,
    None is returned.

    Args:
        pos_sets: Ordered tuple of per-position frozen sets of bases (e.g., {'A', 'G'}).

    Returns:
        Optional[str]: The concatenated IUPAC sequence, or None if a set is not representable.
    """
    try:
        return "".join(IUPAC_LOOKUP[s] for s in pos_sets)
    except KeyError:
        return None  # a set not representable by provided IUPAC map


def compute_pos_sets_from_indices(
    indices: List[int], L: int, pow4: List[int]
) -> Tuple[FrozenSet[str], ...]:
    """
    Compute the set of possible bases at each position across integer-encoded sequences.

    Given integer indices encoding L-length sequences in base-4 (with positional
    weights pow4), extracts the base digit at each position and accumulates the
    union of bases using INT_TO_BASE.

    Args:
        indices: Integers encoding sequences over {A, C, G, T}.
        L: Sequence length (number of positions).
        pow4: Positional weights, where pow4[pos] == 4**pos.

    Returns:
        A tuple of length L where each element is a frozenset of bases observed
        at that position across all indices.
    """
    pos_sets: List[Set[str]] = [set() for _ in range(L)]
    for x in indices:
        for pos in range(L):
            d = (x // pow4[pos]) % 4
            pos_sets[pos].add(INT_TO_BASE[d])
    return tuple(frozenset(s) for s in pos_sets)


def enumerate_product_indices(
    pos_sets: Tuple[FrozenSet[str], ...], pow4: List[int]
) -> List[int]:
    """
    Enumerate linear indices for all sequences in the Cartesian product of per-position base choices.

    Each base is mapped to a digit via BASE_TO_INT (A=0, C=1, G=2, T=3), and indices are computed as
    sum(d_i * pow4[i]) for position i.

    Args:
        pos_sets: Tuple of frozensets with allowed bases at each position.
        pow4: Positional weights, typically 4**i aligned with pos_sets.

    Returns:
        List of integer indices for all combinations. Order reflects set iteration and is not guaranteed.
    """
    # Iterative enumeration to avoid recursion overhead
    idxs = [0]
    for pos, letters in enumerate(pos_sets):
        add = []
        digits = [BASE_TO_INT[b] for b in letters]
        w = pow4[pos]
        for base_idx in idxs:
            for d in digits:
                add.append(base_idx + d * w)
        idxs = add
    return idxs


def is_exact_rectangle(
    pos_sets: Tuple[FrozenSet[str], ...], present: List[bool], pow4: List[int]
) -> bool:
    """
     Check whether a Cartesian product of per-position base choices is fully covered by the presence mask.

     This enumerates linear indices for all sequences in the product using BASE_TO_INT (A=0, C=1, G=2, T=3)
     and positional weights pow4 (typically 4**i). Returns True iff the product is non-empty and every
     corresponding index is marked True in present. By construction, the product contains only intended members.

    Args:
         pos_sets: Tuple of frozensets with allowed bases at each position.
         present: Boolean mask over all 4**L sequences, indexed by sum(d_i * pow4[i]).
         pow4: Positional weights aligned with pos_sets.

     Returns:
         True if all sequences in the product are present; otherwise False.
    """
    idxs = enumerate_product_indices(pos_sets, pow4)
    # Check all in product are present
    if not idxs:
        return False
    for idx in idxs:
        if not present[idx]:
            return False
    return True  # exactness with respect to input set; product contains only members by construction


def close_rectangle(
    pos_sets: Tuple[FrozenSet[str], ...], present: List[bool], L: int, pow4: List[int]
) -> Tuple[FrozenSet[str], ...]:
    """Greedily enlarge a Cartesian product of per-position base choices while preserving exactness.

    This iteratively proposes adding each missing base (A, C, G, T) to each position and accepts an addition only
    if the resulting product remains fully covered in the presence mask, as validated by is_exact_rectangle. The
    process repeats until no further additions are possible, never removing bases.

    Args:
        pos_sets: Tuple of frozensets specifying allowed bases at each position.
        present: Boolean mask over all 4**L sequences, indexed by sum(d_i * pow4[i]).
        L: Number of positions (should equal len(pos_sets)).
        pow4: Positional weights aligned with pos_sets for linear indexing.

    Returns:
        Tuple of frozensets representing the maximal exact rectangle (closure) under the presence constraint.
    """

    pos_sets = [set(s) for s in pos_sets]
    changed = True
    while changed:
        changed = False
        for pos in range(L):
            if len(pos_sets[pos]) == 4:
                continue
            current = pos_sets[pos]
            for b in ("A", "C", "G", "T"):
                if b in current:
                    continue
                trial = tuple(
                    frozenset(s if i != pos else (set(s) | {b}))
                    for i, s in enumerate(pos_sets)
                )
                if is_exact_rectangle(trial, present, pow4):
                    pos_sets[pos].add(b)
                    changed = True
    return tuple(frozenset(s) for s in pos_sets)


# Subset preselection


def preselect_rectangles_by_pairs(
    seqs: Iterable[str], D: int = 2
) -> List[Tuple[str, Set[str]]]:
    """
    Preselect exact degenerate rectangles from equal-length DNA sequences using Hamming-neighbor pairs as seeds.

    Builds a base-4 presence table, seeds from neighbors within distance <= D, verifies all recombinants, greedily closes to a maximal exact rectangle, deduplicates by signature, and also includes singleton rectangles.

    Args:
        seqs (Iterable[str]): Equal-length sequences over A,C,G,T.
        D (int, optional): Maximum Hamming distance for seeding (capped at 2). Defaults to 2.

    Returns:
        List[Tuple[str, Set[str]]]: Tuples of (IUPAC-degenerate code, set of covered input sequences).

    Raises:
        ValueError: If lengths differ or any sequence contains non-ACGT characters.
    """
    # make unique
    seqs = list(set(seqs))
    if not seqs:
        return []
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("All sequences must be the same length")
    if any(ch not in "ACGT" for s in seqs for ch in s):
        raise ValueError("Sequences must contain only A,C,G,T")

    # Precompute powers of 4 (least significant position at index 0)
    pow4 = [1] * L
    for i in range(1, L):
        pow4[i] = pow4[i - 1] * 4

    # Encode sequences and build presence table
    idx_to_seq: Dict[int, str] = {}
    seq_indices: List[int] = []
    max_space = pow4[-1] * 4
    present = [False] * max_space
    for s in seqs:
        x = encode_seq(s, pow4)
        if not present[x]:
            present[x] = True
            idx_to_seq[x] = s
            seq_indices.append(x)

    # Helper to compute minimal rectangle pos_sets from two indices
    def minimal_pos_sets(a: int, b: int) -> Tuple[FrozenSet[str], ...]:
        sets: List[FrozenSet[str]] = []
        for pos in range(L):
            da = (a // pow4[pos]) % 4
            db = (b // pow4[pos]) % 4
            if da == db:
                sets.append(frozenset({INT_TO_BASE[da]}))
            else:
                sets.append(frozenset({INT_TO_BASE[da], INT_TO_BASE[db]}))
        return tuple(sets)

    # Dedup store
    seen: Set[Tuple[FrozenSet[str], ...]] = set()
    out: List[Tuple[str, Set[str]]] = []

    # Include singletons as safe fallbacks
    for x in seq_indices:
        pos_sets = tuple(
            frozenset({INT_TO_BASE[(x // pow4[pos]) % 4]}) for pos in range(L)
        )
        seen.add(pos_sets)
        code = iupac_from_pos_sets(pos_sets)
        if code is not None:
            out.append((code, {idx_to_seq[x]}))

    # Explore neighbor pairs
    for x in seq_indices:
        # Precompute digits of x for speed
        digits_x = [(x // pow4[pos]) % 4 for pos in range(L)]
        for y in neighbors_within_hamming(x, L, pow4, digits_x, D):
            if y <= x or not present[y]:
                continue
            # Build minimal rectangle from the pair
            pos_sets = minimal_pos_sets(x, y)
            # Quick fail: if any position set not representable by provided IUPAC map
            if any(frozenset(s) not in IUPAC_LOOKUP for s in pos_sets):
                continue
            # Check exact rectangle for the 2^d recombinants
            if not is_exact_rectangle(pos_sets, present, pow4):
                continue
            # Greedy closure (try to add third/fourth letters where possible)
            closed = close_rectangle(pos_sets, present, L, pow4)
            if closed in seen:
                continue
            seen.add(closed)
            code = iupac_from_pos_sets(closed)
            if code is None:
                continue
            # Collect covered sequences
            idxs = enumerate_product_indices(closed, pow4)
            covered = {idx_to_seq[i] for i in idxs if present[i]}
            out.append((code, covered))

    return out


## Solving
def minimal_exact_cover(
    seqs: List[str], candidates: List[Tuple[str, Set[str]]], time_limit: float = 30.0
) -> List[str]:
    """
    Solve a minimal exact cover over sequences using a CP-SAT model and degenerate IUPAC codes.

    Each sequence must be covered exactly once by a chosen candidate. The objective is minimized in tie-break order:
    - Minimize number of selected patterns
    - Minimize number of degenerate positions across selected patterns
    - Minimize total allowed bases across all positions

    Args:
        seqs: List of sequences to cover.
        candidates: List of (code, covered_set) where code is the degenerate pattern and covered_set is the set of sequences it covers.
        time_limit: Maximum solver time in seconds.

    Returns:
        List of selected pattern codes forming an exact cover.

    Raises:
        RuntimeError: If no feasible/optimal solution is found within the time limit.
    """
    # Map each sequence to an index
    seq_index = {s: i for i, s in enumerate(seqs)}
    n_seqs = len(seqs)
    n_cands = len(candidates)

    # Build incidence: for each sequence, which candidates cover it
    covers = [[] for _ in range(n_seqs)]
    for cand_idx, (_, covered_set) in enumerate(candidates):
        for s in covered_set:
            covers[seq_index[s]].append(cand_idx)

    # CP-SAT model
    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x_{i}") for i in range(n_cands)]

    # Exact cover constraints: each sequence covered exactly once
    for seq_idx in range(n_seqs):
        model.Add(sum(x[c] for c in covers[seq_idx]) == 1)

    # --- Tie-breaker weights ---
    deg_positions = []
    complexities = []
    for code, _ in candidates:
        deg_positions.append(sum(1 for ch in code if len(IUPAC_CODES[ch]) > 1))
        complexities.append(sum(len(IUPAC_CODES[ch]) for ch in code))

    # Weight scales: BIG1 >> BIG2 >> BIG3
    BIG1 = 10**6  # dominates everything else
    BIG2 = 10**3  # dominates complexity term

    model.Minimize(
        sum(x[i] * BIG1 for i in range(n_cands))  # primary: #patterns
        + sum(
            x[i] * BIG2 * deg_positions[i] for i in range(n_cands)
        )  # secondary: degenerate positions
        + sum(x[i] * complexities[i] for i in range(n_cands))  # tertiary: total bases
    )

    # Solve
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8

    status = solver.Solve(model)

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError("No solution found")

    # Extract chosen patterns
    chosen_codes = [candidates[i][0] for i in range(n_cands) if solver.Value(x[i]) == 1]
    return chosen_codes

def seqs_to_degenerates(seqs: Iterable[str], D: int = 2, time_limit: float = 30.0) -> List[str]:
    """Takes a list of equal length sequences and finds a minimal set of degenerate codes that represents them.
    Args:
        seqs: A list of sequences where each sequence is a string of equal length. The sequences should only contain valid nucleotide characters (e.g., A, T, C, G).
    Returns:
        a list of degenerate sequences representing the input list
    """
    seqs = list(set(seqs))
    if not seqs:
        return []
    elif len(seqs) == 1:
        return seqs
    try:
        candidates = preselect_rectangles_by_pairs(seqs=seqs, D=D)
        return minimal_exact_cover(seqs=seqs, candidates=candidates, time_limit=time_limit)
    except Exception as e:
        logging.exception(f"An error occurred: {e}")
        raise e
    