"""Pascal's triangle mod n as a fractal, and what its geometry says.

Plotted, Pascal's triangle mod a prime ``p`` is a Sierpinski gasket with scaling
ratio ``p``.  That is not decoration -- it *is* Lucas' theorem: the entry at
``(i, j)`` is non-zero mod ``p`` exactly when every base-``p`` digit of ``j`` is
at most the corresponding digit of ``i``, which is precisely the self-similar
construction rule of the gasket.

For composite ``n = pq`` the picture is two gaskets superimposed, at scales ``p``
and ``q``.  The support is the *union* of the two supports (an entry vanishes mod
``n`` only when it vanishes mod both), so the zeros of the composite picture are
the places where the two hole systems coincide.

This module makes that geometry computable:

``row_nonzero_count``   Lucas' count of surviving entries in one row.
``support_count``       exact box-count over the first ``N`` rows, in ``O(log N)``.
``fractal_dimension``   ``log(p(p+1)/2) / log p``, the gasket's box dimension.
``first_zero_row``      the row where the two hole systems first coincide.

The geometric reading also explains the project's central barrier in one
sentence: **you cannot see a scale-``p`` fractal by sampling below scale ``p``.**
Theorem 7's aliasing bound is that statement in arithmetic dress.
"""

from __future__ import annotations

from math import log, prod

__all__ = [
    "digits",
    "column_divides",
    "row_nonzero_count",
    "row_has_zero",
    "support_count",
    "fractal_dimension",
    "first_zero_row",
    "first_zero_row_predicted",
    "render",
]


def digits(n: int, base: int) -> list[int]:
    """Base-``base`` digits of ``n``, least significant first."""
    if n == 0:
        return [0]
    out = []
    while n:
        out.append(n % base)
        n //= base
    return out


def row_nonzero_count(i: int, p: int) -> int:
    """How many entries of row ``i`` survive mod ``p``.

    Lucas: the survivors are the ``j`` whose base-``p`` digits are dominated by
    those of ``i``, so the count is ``prod (d + 1)`` over the digits of ``i``.
    """
    return prod(d + 1 for d in digits(i, p))


def row_has_zero(i: int, p: int) -> bool:
    """Does row ``i`` contain an entry divisible by ``p``?"""
    return row_nonzero_count(i, p) < i + 1


def support_count(rows: int, p: int) -> int:
    """Exact number of entries not divisible by ``p`` in rows ``0..rows-1``.

    A digit dynamic program, so this is ``O(log_p rows)`` -- the box-count of the
    gasket without drawing it.  For ``rows = p**k`` it returns
    ``(p(p+1)/2)**k``, the classical self-similarity relation.
    """
    if rows <= 0:
        return 0
    nd = digits(rows, p)
    block = p * (p + 1) // 2
    total = 0
    # digits above position t are fixed to those of `rows`; digit t is smaller;
    # everything below t is free.
    for t in range(len(nd) - 1, -1, -1):
        above = prod(nd[s] + 1 for s in range(t + 1, len(nd)))
        here = nd[t] * (nd[t] + 1) // 2
        total += above * here * block**t
    return total


def fractal_dimension(p: int) -> float:
    """Box dimension of the mod-``p`` gasket, ``log(p(p+1)/2) / log p``."""
    return log(p * (p + 1) / 2) / log(p)


def first_zero_row_predicted(p: int, q: int, cap: int | None = None) -> int | None:
    """Least ``i`` whose row contains a zero mod ``p*q``, from the digit criterion.

    An entry vanishes mod ``pq`` only when it vanishes mod both, so the row must
    carry a hole in *both* gaskets.  Verified against direct computation.
    """
    cap = cap or 6 * max(p, q) + 60
    for i in range(2, cap):
        if row_has_zero(i, p) and row_has_zero(i, q):
            return i
    return None


def first_zero_row(n: int, cap: int | None = None) -> int | None:
    """Least row of Pascal's triangle containing an entry ``== 0 mod n``.

    Computed directly, row by row.  Costs ``O(cap**2)`` operations -- reading the
    two-dimensional picture costs *area*, which is why this is not a shortcut.
    """
    cap = cap or 4 * n
    row = [1]
    for i in range(1, cap):
        row = [1] + [(row[j - 1] + row[j]) % n for j in range(1, i)] + [1]
        if any(v == 0 for v in row[1:-1]):
            return i
    return None


def render(rows: int, n: int, legend: tuple[int, int] | None = None) -> list[str]:
    """ASCII picture of the triangle mod ``n``.

    With ``legend = (p, q)`` the two gaskets are distinguished: ``p`` marks an
    entry killed by ``p`` alone, ``q`` by ``q`` alone, ``0`` by both (a genuine
    zero mod ``n``), ``#`` a survivor.
    """
    out = []
    if legend is None:
        row = [1]
        for i in range(rows):
            if i:
                row = [1] + [(row[j - 1] + row[j]) % n for j in range(1, i)] + [1]
            out.append(" " * (rows - i - 1) + "".join("#" if v else "." for v in row))
        return out
    p, q = legend
    rp, rq = [1], [1]
    for i in range(rows):
        if i:
            rp = [1] + [(rp[j - 1] + rp[j]) % p for j in range(1, i)] + [1]
            rq = [1] + [(rq[j - 1] + rq[j]) % q for j in range(1, i)] + [1]
        line = " " * (rows - i - 1)
        for j in range(i + 1):
            zp, zq = rp[j] == 0, rq[j] == 0
            line += "0" if (zp and zq) else ("p" if zp else ("q" if zq else "#"))
        out.append(line)
    return out


def column_divides(i: int, c: int, p: int) -> bool:
    """Does ``p`` divide ``C(i, c)``?

    Lucas: ``p`` divides ``C(i, c)`` exactly when some base-``p`` digit of ``c``
    exceeds the corresponding digit of ``i``.  Reading *down* column ``c`` this
    is a condition on ``i`` modulo a power of ``p``, so the column's zero pattern
    is periodic -- the sideways view of the gasket.
    """
    dc, di = digits(c, p), digits(i, p)
    for t, x in enumerate(dc):
        y = di[t] if t < len(di) else 0
        if x > y:
            return True
    return False
