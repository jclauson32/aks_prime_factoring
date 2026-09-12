"""Round 17: Pascal's triangle on the collision mechanism.

Iterate a column of the triangle as a map on ``Z/N``:

    x  ->  C(x, k) mod N.

For ``k = 2`` this *is* Pollard's rho, up to an affine change of variable:
with ``x = 2y + 1/2``,

    C(x, 2)  =  2 (y^2 - 5/16) + 1/2,

so iterating ``C(x, 2)`` is iterating ``y^2 - 5/16`` (``conjugate_check``).
For larger ``k`` the map is a degree-``k`` polynomial, and how fast its orbits
collide mod ``p`` is governed by one number, the fibre statistic

    kappa = (1/p) sum_v m_v (m_v - 1),        m_v = #{x : C(x, k) = v},

with expected rho length ``sqrt(pi p / (2 kappa))``.  Pascal's reflection
``C(k-1-x, k) = (-1)^k C(x, k)`` decides it: for even ``k >= 4`` the map factors
through ``(x - (k-1)/2)^2``, every fibre is closed under the reflection, and
``kappa = 2``; for ``k = 2`` and odd ``k`` it is ``1``, as for a random map.

Concretely, with ``u = x - 3/2``, ``C(x, 4) = ((u^2 - 5/4)^2 - 1) / 24``: two
squarings per step buy a factor ``sqrt 2`` in steps, which is what the second
squaring costs.
"""

from __future__ import annotations

from collections import Counter
from math import factorial, gcd


def binom_poly(k: int, n: int):
    """The map ``x -> C(x, k) mod n`` (needs ``gcd(k!, n) = 1``)."""
    inv = pow(factorial(k), -1, n)

    def f(x: int) -> int:
        r = inv
        for j in range(k):
            r = r * (x - j) % n
        return r

    return f


def conjugate_check(p: int, limit: int = 5000) -> bool:
    """Check ``C(2y + 1/2, 2) = 2(y^2 - 5/16) + 1/2`` on ``F_p``."""
    half, sixteenth = pow(2, -1, p), pow(16, -1, p)
    for y in range(min(p, limit)):
        x = (2 * y + half) % p
        if x * (x - 1) * half % p != (2 * (y * y - 5 * sixteenth) + half) % p:
            return False
    return True


def fibre_kappa(k: int, p: int) -> float:
    """``(1/p) sum m_v (m_v - 1)`` for ``x -> C(x, k)`` on ``F_p``."""
    counts = Counter(map(binom_poly(k, p), range(p)))
    return sum(m * (m - 1) for m in counts.values()) / p


def predicted_kappa(k: int) -> int:
    """``2`` for even ``k >= 4`` (reflection-symmetric fibres), else ``1``."""
    return 2 if k >= 4 and k % 2 == 0 else 1


def rho_length(f, x0: int) -> int:
    """Tail plus cycle length of the orbit of ``x0`` under ``f`` (small moduli)."""
    seen = {}
    x = x0
    while x not in seen:
        seen[x] = len(seen)
        x = f(x)
    return len(seen)


def degenerate_start(k: int, x0: int) -> bool:
    """Whether ``x0`` is one of the starts the map cannot walk away from.

    Over the integers ``x -> C(x, k)`` has two fixed points, ``0`` and
    ``k + 1`` (``C(k+1, k) = k+1``), and every ``x0 <= k`` lands on ``0`` or
    ``1`` within two steps and then stays at ``0``.  A walk from such a start
    stops moving, so ``x - y`` is ``0``, every gcd is ``n``, and the search
    reports failure after a handful of evaluations.  For ``k = 1`` the map is
    the identity and no start walks at all.
    """
    return k <= 1 or x0 <= k + 1


def pascal_rho(n: int, k: int = 2, x0: int | None = None, max_steps: int = 1 << 26):
    """Brent's cycle-finding on ``x -> C(x, k) mod n``, batching gcds.

    Returns ``(factor, steps)`` where ``steps`` counts every evaluation of the
    map, including Brent's advance loop; ``factor`` is ``None`` if the walk
    closed up mod ``n`` (retry with another ``x0``) or ran out of steps.

    ``x0`` defaults to ``k + 2``, the smallest start that is not degenerate for
    the given ``k`` (see ``degenerate_start``).
    """
    if x0 is None:
        x0 = k + 2
    if n % 2 == 0:
        return 2, 0
    g = gcd(factorial(k), n)
    if 1 < g < n:
        return g, 0
    f = binom_poly(k, n)
    y, r, q, steps = x0 % n, 1, 1, 0
    batch = 64
    while steps < max_steps:
        x = y
        for _ in range(r):
            y = f(y)
        steps += r
        done = 0
        while done < r:
            ys = y
            for _ in range(min(batch, r - done)):
                y = f(y)
                q = q * (x - y) % n
            steps += min(batch, r - done)
            done += min(batch, r - done)
            g = gcd(q, n)
            if g > 1:
                if g == n:
                    # back up and step singly
                    while True:
                        ys = f(ys)
                        steps += 1
                        g = gcd(x - ys, n)
                        if g > 1:
                            break
                return (g if g < n else None), steps
        r *= 2
    return None, steps
