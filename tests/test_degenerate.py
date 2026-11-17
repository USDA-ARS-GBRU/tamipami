import random

import pytest

from tamipami import degenerate
from tamipami.candidates import IUPAC_TO_BASES, expand_iupac

DEGENERATE_SYMBOLS = sorted(IUPAC_TO_BASES)


def _random_degenerate_codes(k, n, seed):
    """n distinct random IUPAC codes of length k with at least one degenerate position."""
    rng = random.Random(seed)
    codes = set()
    while len(codes) < n:
        code = "".join(rng.choice("ACGT" + "".join(DEGENERATE_SYMBOLS)) for _ in range(k))
        if any(len(IUPAC_TO_BASES[c]) > 1 for c in code):
            codes.add(code)
    return sorted(codes)


@pytest.mark.parametrize("k", [3, 4, 5])
def test_single_code_expansion_is_recovered_as_one_code(k):
    """The full expansion of one degenerate code is covered by exactly that code."""
    for code in _random_degenerate_codes(k, 25, seed=k):
        result = degenerate.solve_degenerates(expand_iupac(code))
        assert result.codes == [code]
        assert result.proven_optimal
