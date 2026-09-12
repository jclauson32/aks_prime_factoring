"""Round 32: Shanks' square forms factorisation -- the sign mechanism without smoothness.

SQUFOF walks the cycle of reduced binary quadratic forms of discriminant ``4N``
(equivalently the continued fraction of ``sqrt(N)``) until it meets a *square*
form.  Taking its square root and walking back gives an ambiguous form, whose
outer coefficient is a factor of ``N``.  It is a congruence of squares -- the
same mechanism as the quadratic sieve and as round 15's central column -- but
it finds the congruence by walking, not by collecting smooth relations, and it
costs ``O(N^(1/4))``.

That is the point of having it here: the sign mechanism *without* smoothness
sits at ``N^(1/4)``, exactly where the central column sits.  What buys the
sieves their sub-exponential time is the smooth numbers, not the squares.
"""

from __future__ import annotations

from math import gcd, isqrt

MULTIPLIERS = (1, 3, 5, 7, 11, 3 * 5, 3 * 7, 3 * 11, 5 * 7, 5 * 11, 7 * 11,
               3 * 5 * 7, 3 * 5 * 11, 3 * 7 * 11, 5 * 7 * 11, 3 * 5 * 7 * 11)


def _squfof_one(n: int, max_steps: int):
    """One pass of SQUFOF on ``n`` (already multiplied by a multiplier)."""
    s = isqrt(n)
    if s * s == n:
        return s
    q0, p0, q1 = 1, s, n - s * s
    if q1 == 0:
        return s
    queue = []
    b = 0
    for i in range(1, max_steps):
        b = (s + p0) // q1
        p1 = b * q1 - p0
        q2 = q0 + b * (p0 - p1)
        if i % 2 == 0:                      # even index: q1 may be a square
            r = isqrt(q1)
            if r * r == q1 and r not in queue:
                # walk back from the square form
                qq0 = r
                pp = p0 + r * ((s - p0) // r)
                qq1 = (n - pp * pp) // qq0
                while True:
                    bb = (s + pp) // qq1
                    pp2 = bb * qq1 - pp
                    if pp2 == pp:
                        break
                    qq2 = qq0 + bb * (pp - pp2)
                    qq0, qq1, pp = qq1, qq2, pp2
                g = gcd(n, pp)
                if 1 < g < n:
                    return g
        if i % 2 == 0:
            queue.append(isqrt(q1) if isqrt(q1) ** 2 == q1 else 0)
        q0, p0, q1 = q1, p1, q2
        if q1 == 0:
            break
    return None


def squfof(n: int, max_steps: int | None = None, stats: dict | None = None):
    """A non-trivial factor of the odd composite ``n``, or ``None``.

    Tries the standard small multipliers; each pass walks the principal cycle
    for ``O(N^(1/4))`` steps looking for a square form.
    """
    if n % 2 == 0:
        return 2
    r = isqrt(n)
    if r * r == n:
        return r
    limit = max_steps or (2 * isqrt(2 * isqrt(n)) + 40)
    for k, mult in enumerate(MULTIPLIERS, start=1):
        if mult * n > 1 << 512:
            break
        g = _squfof_one(mult * n, limit)
        if g:
            g = gcd(g, n)
            if 1 < g < n:
                if stats is not None:
                    stats.update(multiplier=mult, multipliers_tried=k)
                return g
    if stats is not None:
        stats.update(multipliers_tried=len(MULTIPLIERS))
    return None
