"""Round 18: factoring by counting lattice points under a hyperbola.

The number of divisors of ``N`` up to ``X`` is a difference of two lattice-point
counts:

    #{d | N : d <= X}  =  sum_{y <= X} floor(N/y)  -  sum_{y <= X} floor((N-1)/y),

and it is monotone in ``X`` -- the binary-searchable predicate asked for in the
"is there a binary search?" question, in a form with no factorials in it.

Each count is computed exactly by walking the convex hull of the lattice points
*above* the hyperbola ``xy = N``.  That region is convex, so no lattice point lies
between its hull and the curve, and the count in row ``y`` is read off the hull
edge crossing that row.  The hull has ``Theta(N^(1/3) log N)`` vertices between
rows ``N^(1/3)`` and ``N^(1/2)`` (curvature: each dyadic block of the arc has about
``N^(1/3)``), and the Stern-Brocot stack finds each one in amortised ``O(1)``
steps -- the classical ``O(N^(1/3) log N)`` method for the divisor problem
(Vinogradov; made practical by Sladkey).

Walking the hulls for ``N`` and ``N - 1`` side by side, the cumulative row counts
agree until the first divisor, so one pair of walks locates ``p``: factoring in
``O~(N^(1/3))``, by a route that has neither a group, a factorial, nor a square
root in it.  The edges are the exact linear pieces of ``floor(N/y)``; a
faster method would need exact *curved* pieces (``hyperbola_piece_exponent``),
and already for parabolas the full-period floor sum computes class numbers
(``quadratic_floor_sum``).
"""

from __future__ import annotations

import bisect
from math import gcd, isqrt


def icbrt(n: int) -> int:
    """Integer cube root."""
    r = int(round(n ** (1.0 / 3.0))) if n < 1 << 900 else 1 << (n.bit_length() // 3)
    while r * r * r > n:
        r -= 1
    while (r + 1) ** 3 <= n:
        r += 1
    return r


def _edge_rows(x: int, dx: int, dy: int) -> int:
    """Lattice points strictly left of the edge ``(x, y) -> (x+dx, y-dy)`` in its
    ``dy`` rows ``y-dy .. y-1``: ``sum_{j=1}^{dy} (x - 1 + ceil(j dx / dy))``."""
    return dy * x + (dy + 1) * (dx - 1) // 2 if dy else 0


def row_sum(n: int, top: int, stats: dict | None = None):
    """``sum_{y=1}^{top} floor(n/y)`` exactly, for ``1 <= top <= isqrt(n)``.

    Returns ``(total, checkpoints)``; each checkpoint ``(r, t)`` records
    ``t = sum_{y=r}^{top} floor(n/y)`` at a hull vertex in row ``r``.
    """
    if not 1 <= top <= isqrt(n):
        raise ValueError("need 1 <= top <= isqrt(n)")

    def outside(x, y):
        return x * y > n

    y = top
    x = n // top + 1                      # leftmost lattice point above the curve
    total = n // top
    checkpoints = [(y, total)]
    stack = [(1, 0), (0, 1)]              # Stern-Brocot neighbours, steepest on top
    low = icbrt(n)
    steps = 0
    while y > low:
        dx1, dy1 = stack.pop()
        while y - dy1 >= 1 and outside(x + dx1, y - dy1):
            total += _edge_rows(x, dx1, dy1)
            x += dx1
            y -= dy1
            steps += 1
            if dy1:
                checkpoints.append((y, total))
        if y <= low:
            break
        dx2, dy2 = dx1, dy1
        while stack:
            dx1, dy1 = stack[-1]
            if outside(x + dx1, y - dy1):
                break
            dx2, dy2 = dx1, dy1
            stack.pop()
        if not stack:
            break
        while True:                       # refine between feasible (1) and infeasible (2)
            mx, my = dx1 + dx2, dy1 + dy2
            steps += 1
            if outside(x + mx, y - my):
                stack.append((mx, my))
                dx1, dy1 = mx, my
            elif (x + mx) * (x + mx) * dy1 >= n * dx1:
                break                     # (dx1, dy1) is at least as steep as the curve there
            else:
                dx2, dy2 = mx, my
    for r in range(y - 1, 0, -1):
        total += n // r
        steps += 1
    if stats is not None:
        stats["steps"] = steps
        stats["vertices"] = len(checkpoints)
    return total, checkpoints


def row_sum_naive(n: int, top: int) -> int:
    """Reference: ``sum_{y<=top} floor(n/y)`` by quotient blocks, ``O(sqrt n)``."""
    s, y = 0, 1
    while y <= top:
        q = n // y
        y2 = min(top, n // q)
        s += q * (y2 - y + 1)
        y = y2 + 1
    return s


def divisors_up_to(n: int, x: int) -> int:
    """``#{d | n : d <= x}`` for ``x <= isqrt(n - 1)``, from two hyperbola counts."""
    return row_sum(n, x)[0] - row_sum(n - 1, x)[0]


def hyperbola_factor(n: int, stats: dict | None = None):
    """Largest divisor of ``n`` in ``(cbrt n, sqrt n]``, or ``None``, by walking the
    hulls of ``xy > n`` and ``xy > n-1`` and finding the first row where the
    cumulative counts differ.  Divisors below ``cbrt n`` are left to trial
    division (``O(n^(1/3))``, the same budget)."""
    top = isqrt(n - 1)
    sa, sb = {}, {}
    _, a = row_sum(n, top, sa)
    _, b = row_sum(n - 1, top, sb)
    if stats is not None:
        stats["steps"] = sa["steps"] + sb["steps"]
        stats["vertices"] = sa["vertices"] + sb["vertices"]
    rows_b = [r for r, _ in b][::-1]
    tot_b = [t for _, t in b][::-1]

    def cum_b(r):
        i = bisect.bisect_left(rows_b, r)
        if rows_b[i] == r:
            return tot_b[i]
        return tot_b[i] + sum((n - 1) // y for y in range(r, rows_b[i]))

    prev = top + 1
    for r, t in a:
        if t != cum_b(r):
            for y in range(prev - 1, r - 1, -1):
                if n % y == 0:
                    return y
        prev = r
    return None


def hyperbola_piece_exponent(degree: int) -> float:
    """Pieces needed if ``floor(N/y)`` were summed with exact degree-``d`` arcs.

    On ``[x, 2x]`` a degree-``d`` Taylor piece of length ``h`` misses ``N/y`` by
    ``~ N h^(d+1) / x^(d+2)``; asking for ``O(1)`` stray lattice points per piece
    gives ``h ~ x N^(-1/(d+2))``, i.e. ``N^(1/(d+2))`` pieces per dyadic block.
    ``d = 1`` is the hull walk, ``N^(1/3)``.
    """
    return 1.0 / (degree + 2)


def quadratic_floor_sum(p: int) -> int:
    """``sum_{k=0}^{p-1} floor(k^2 / p)``.  For primes ``p = 3 (mod 4)``, ``p > 3``,
    this equals ``(p-1)(2p-1)/6 - (p-1-2h)/2`` with ``h = h(-p)`` the class number
    (Dirichlet: ``h(-p) = -(1/p) sum a (a/p)``), so an exact fast quadratic floor
    sum would compute class numbers."""
    return sum(k * k // p for k in range(p))


def class_number(d: int) -> int:
    """Class number of the negative discriminant ``d`` by counting reduced forms."""
    if d >= 0 or d % 4 not in (0, 1):
        raise ValueError("need a negative discriminant")
    h, a = 0, 1
    while 3 * a * a <= -d:
        for b in range(-a + 1, a + 1):
            if (b * b - d) % (4 * a):
                continue
            c = (b * b - d) // (4 * a)
            if c < a or (b < 0 and a == c):
                continue
            if gcd(gcd(a, abs(b)), c) == 1:
                h += 1
        a += 1
    return h
