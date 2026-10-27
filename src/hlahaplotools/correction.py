"""Multiple-testing correction."""

from __future__ import annotations

from collections.abc import Sequence


def adjust_fdr(p_values: Sequence[float | None]) -> list[float | None]:
    """Benjamini–Hochberg adjusted P values (q-values), in the input order.

    Missing values (None) are skipped and stay None; the number of tests is the number
    of P values given.
    """
    idx = [i for i, p in enumerate(p_values) if p is not None]
    m = len(idx)
    out: list[float | None] = [None] * len(p_values)
    running = 1.0
    for rank, i in enumerate(sorted(idx, key=lambda i: p_values[i], reverse=True)):
        q = p_values[i] * m / (m - rank)
        running = min(running, q)
        out[i] = running
    return out
