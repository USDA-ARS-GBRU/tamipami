#!/usr/bin/env python
"""Independent global-optimum oracle for minimal exact degenerate covers.

The oracle shares no code with tamipami's candidate generation. It enumerates every
IUPAC code of length k over the full alphabet (15^k codes), keeps the codes whose
expansion is a subset of the input set D, and solves the three-level lexicographic
exact-cover problem (patterns, degenerate positions, total bases) to proven
optimality with CP-SAT.

`tamipami_solve` runs tamipami.degenerate.seqs_to_degenerates on the same input,
checks that the result is an exact cover, and records whether every CP-SAT stage
was proven optimal and whether candidate enumeration was limited.
"""
import itertools
import logging
from typing import List, Optional, Tuple

from ortools.sat.python import cp_model

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "M": "AC", "R": "AG", "W": "AT", "S": "CG", "Y": "CT", "K": "GT",
    "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "N": "ACGT",
}

Score = Tuple[int, int, int]


def expand(code: str) -> List[str]:
    return ["".join(t) for t in itertools.product(*(IUPAC[c] for c in code))]


def score(codes: List[str]) -> Score:
    """(#patterns, #degenerate positions, #bases), the lexicographic objective."""
    return (
        len(codes),
        sum(len(IUPAC[c]) > 1 for code in codes for c in code),
        sum(len(IUPAC[c]) for code in codes for c in code),
    )


def exact_codes(D: List[str]) -> List[Tuple[str, List[str]]]:
    """Every IUPAC code of length k whose expansion is contained in D."""
    d_set, k = set(D), len(D[0])
    out = []
    for code in itertools.product(IUPAC, repeat=k):
        exp = []
        for t in itertools.product(*(IUPAC[c] for c in code)):
            s = "".join(t)
            if s not in d_set:
                break
            exp.append(s)
        else:
            out.append(("".join(code), exp))
    return out


def oracle(D: List[str], time_limit: float = 300.0, workers: int = 8) -> Tuple[Optional[Score], int]:
    """
    Return (globally optimal score, number of exact codes), or (None, n) if any
    objective level is not proven optimal within `time_limit` seconds.
    """
    cands = exact_codes(D)
    idx = {s: i for i, s in enumerate(D)}
    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x{j}") for j in range(len(cands))]
    covers = [[] for _ in D]
    for j, (_, exp) in enumerate(cands):
        for s in exp:
            covers[idx[s]].append(j)
    for c in covers:
        model.AddExactlyOne(x[j] for j in c)
    objectives = [
        cp_model.LinearExpr.Sum(x),
        cp_model.LinearExpr.WeightedSum(x, [sum(len(IUPAC[ch]) > 1 for ch in c) for c, _ in cands]),
        cp_model.LinearExpr.WeightedSum(x, [sum(len(IUPAC[ch]) for ch in c) for c, _ in cands]),
    ]
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    for obj in objectives:
        model.Minimize(obj)
        if solver.Solve(model) != cp_model.OPTIMAL:
            return None, len(cands)
        model.Add(obj == int(round(solver.ObjectiveValue())))
        model.ClearHints()
        for v in x:
            model.AddHint(v, solver.Value(v))
    return score([cands[j][0] for j in range(len(cands)) if solver.Value(x[j])]), len(cands)


class _StatusLog(logging.Handler):
    """Collects CP-SAT stage statuses and candidate-limit messages from tamipami logs."""

    def __init__(self):
        super().__init__()
        self.statuses: List[str] = []
        self.limited = False

    def emit(self, record):
        msg = record.getMessage()
        if "status=" in msg:
            self.statuses.append(msg.split("status=")[1].split()[0].replace("CpSolverStatus.", ""))
        if "Skipping positional" in msg or "Kept the" in msg:
            self.limited = True


def tamipami_solve(D: List[str]) -> Tuple[Score, bool, bool]:
    """Return (score, proven optimal at every stage, candidate enumeration limited)."""
    from tamipami import degenerate

    handler = _StatusLog()
    loggers = [logging.getLogger(n) for n in ("tamipami.degenerate", "tamipami.candidates")]
    for lg in loggers:
        lg.addHandler(handler)
        lg.setLevel(logging.INFO)
    try:
        codes = degenerate.seqs_to_degenerates(D)
    finally:
        for lg in loggers:
            lg.removeHandler(handler)
    covered = sorted(s for c in codes for s in expand(c))
    if covered != sorted(D):
        raise AssertionError(f"tamipami result is not an exact cover: {codes}")
    proven = not handler.limited and all(st == "OPTIMAL" for st in handler.statuses)
    return score(codes), proven, handler.limited
