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

No binary search is needed.  Every divisor point ``(N/d, d)`` lies on the
strictly convex curve ``xy = N``, so it is a *vertex* of the hull of
``{xy > N - 1}``, and one walk of that hull finds it: factoring in
``O~(N^(1/3))``, by a route that has neither a group, a factorial, nor a square
root in it.  The two hull edges at the divisor are Farey neighbours bracketing
``p/q``, and each is a Lehman certificate (``divisor_vertex_edges``).  The edges are the exact linear pieces of ``floor(N/y)``; a
faster method would need exact *curved* pieces (``hyperbola_piece_exponent``),
and already for parabolas the full-period floor sum computes class numbers
(``quadratic_floor_sum``).
"""

from __future__ import annotations

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


def hull_walk(n: int, top: int, stats: dict | None = None):
    """Yield the lattice points ``(x, y, (dx, dy))`` along the lower-left hull of
    ``{xy > n}``, from row ``top`` down to row ``cbrt(n)``; ``(dx, dy)`` is the
    primitive step that arrived there (``None`` at the start)."""

    def outside(x, y):
        return x * y > n

    y = top
    x = n // top + 1
    yield x, y, None
    stack = [(1, 0), (0, 1)]
    low = icbrt(n)
    steps = 0
    while y > low:
        dx1, dy1 = stack.pop()
        while y - dy1 >= 1 and outside(x + dx1, y - dy1):
            x += dx1
            y -= dy1
            steps += 1
            yield x, y, (dx1, dy1)
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
        while True:
            mx, my = dx1 + dx2, dy1 + dy2
            steps += 1
            if outside(x + mx, y - my):
                stack.append((mx, my))
                dx1, dy1 = mx, my
            elif (x + mx) * (x + mx) * dy1 >= n * dx1:
                break
            else:
                dx2, dy2 = mx, my
        if stats is not None:
            stats["steps"] = steps
    if stats is not None:
        stats["steps"] = steps


def hyperbola_factor(n: int, stats: dict | None = None):
    """Largest divisor of ``n`` in ``[cbrt n, sqrt n)``, or ``None``: one hull walk.

    ``Outside(n-1) = Outside(n) + {(x, y) : xy = n}``, and a lattice point on the
    strictly convex curve ``xy = n`` is an extreme point of ``{xy >= n}``, so
    every divisor point ``(n/d, d)`` is a *vertex* of the hull of ``{xy > n-1}``.
    Walk that hull and test ``x * y == n`` at each point.  Divisors below
    ``cbrt n`` are left to trial division (``O(n^(1/3))``, the same budget).
    """
    st = {}
    for x, y, _ in hull_walk(n - 1, isqrt(n - 1), st):
        if x * y == n:
            if stats is not None:
                stats["steps"] = st.get("steps", 0)
            return y
    if stats is not None:
        stats["steps"] = st.get("steps", 0)
    return None


def divisor_vertex_edges(n: int):
    """The two hull edges of ``{xy > n-1}`` at the divisor vertex ``(q, p)``.

    Each primitive edge ``(dx, dy)`` gives Lehman's identity
    ``(dy q + dx p)^2 - 4 (dx dy) n = (dy q - dx p)^2``; returns
    ``(p, q, [(dx, dy, k, gap), ...])`` with ``k = dx dy``, ``gap = |dy q - dx p|``.
    """
    pts = list(hull_walk(n - 1, isqrt(n - 1)))
    for i, (x, y, e) in enumerate(pts):
        if x * y == n:
            edges = [e]
            if i + 1 < len(pts):
                edges.append(pts[i + 1][2])
            out = [(dx, dy, dx * dy, abs(dy * x - dx * y)) for dx, dy in (e for e in edges if e)]
            return y, x, out
    return None


def vertex_cones(n: int, top: int, stop_row: int):
    """Hull vertices of ``{xy > n}`` with the angle of their normal cones.

    Walks from row ``top`` down to ``stop_row``; returns ``(x, y, angle)`` for
    every point where the edge direction turns, ``angle`` being the turn in
    radians -- the range of slopes whose supporting line touches that vertex.
    """
    from math import atan2

    pts = []
    for pt in hull_walk(n, top):
        pts.append(pt)
        if pt[1] < stop_row:
            break
    out = []
    for i in range(1, len(pts) - 1):
        e_in, e_out = pts[i][2], pts[i + 1][2]
        if e_in and e_out and e_in != e_out:
            out.append((pts[i][0], pts[i][1],
                        abs(atan2(e_out[1], e_out[0]) - atan2(e_in[1], e_in[0]))))
    return out


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
