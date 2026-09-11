"""The central column of Pascal's triangle computes Legendre symbols.

    sum_k C(2k, k) x^k  =  (1 - 4x)^(-1/2)

Modulo a prime ``p``, for ``0 <= k <= (p-1)/2`` one has
``C(2k,k) = C((p-1)/2, k) (-4)^k``, and for ``(p+1)/2 <= k <= p-1`` Lucas makes
``C(2k,k) = 0``.  So the **truncated** column sum

    P_T(x) = sum_{k <= T} C(2k, k) x^k   ==   (1 - 4x)^((p-1)/2)   (mod p)

for *every* ``T`` in the window ``[(p-1)/2, p-1]`` -- and by Euler's criterion
``(1-4a)^((p-1)/2)`` is the Legendre symbol ``((1-4a)/p)``.

That gives a factoring method using a mechanism nothing else in this repository
uses: the **sign** of a square root.  ``P_T(a) mod N`` is ``+-1`` modulo ``p`` and
something else modulo ``q``, so ``gcd(P_T(a) -+ 1, N)`` splits ``N``.  Three
properties make it worth having:

* the window ``[(p-1)/2, p-1]`` spans a factor of two, so a geometric search over
  ``T = 2^j`` is guaranteed to land in it -- no scanning;
* it fires at ``T ~ p/2``, *below* the threshold where Strassen's factorial test
  (Theorem 20) sees anything at all, since ``gcd(T!, N) = 1`` for ``T < p``;
* ``P_T`` is a hypergeometric partial sum, so it costs ``O~(sqrt(T))`` by the
  same baby-step/giant-step as the factorial, run on 2x2 polynomial matrices.

Net cost ``O~(sqrt(p))``: Strassen's exponent, reached through a different door.
"""

from __future__ import annotations

from math import gcd, isqrt

from .fast import factorial_mod, multipoint_eval, poly_mul

__all__ = ["jacobi", "central_sum_naive", "central_sum", "central_split"]


def jacobi(a: int, n: int) -> int:
    """The Jacobi symbol ``(a/n)`` for odd ``n > 0``."""
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def central_sum_naive(t: int, a: int, n: int) -> int:
    """``sum_{k<=t} C(2k,k) a^k mod n`` in ``O(t)``, via the term ratio.

    ``C(2k+2,k+1)/C(2k,k) = 2(2k+1)/(k+1)``; the division needs ``k+1`` coprime
    to ``n``, i.e. ``t < spf(n)``.
    """
    s = term = 1 % n
    for k in range(t):
        term = term * (2 * (2 * k + 1) * a) % n * pow(k + 1, -1, n) % n
        s = (s + term) % n
    return s


def _padd(u, v, n):
    out = [0] * max(len(u), len(v))
    for i, c in enumerate(u):
        out[i] = c
    for i, c in enumerate(v):
        out[i] = (out[i] + c) % n
    while out and out[-1] == 0:
        out.pop()
    return out


def _mmul(a, b, n):
    """2x2 matrices of polynomials over Z/n."""
    return [
        [_padd(poly_mul(a[0][0], b[0][0], n), poly_mul(a[0][1], b[1][0], n), n),
         _padd(poly_mul(a[0][0], b[0][1], n), poly_mul(a[0][1], b[1][1], n), n)],
        [_padd(poly_mul(a[1][0], b[0][0], n), poly_mul(a[1][1], b[1][0], n), n),
         _padd(poly_mul(a[1][0], b[0][1], n), poly_mul(a[1][1], b[1][1], n), n)],
    ]


def _nmul(a, b, n):
    """2x2 numeric matrices over Z/n."""
    return [
        [(a[0][0] * b[0][0] + a[0][1] * b[1][0]) % n,
         (a[0][0] * b[0][1] + a[0][1] * b[1][1]) % n],
        [(a[1][0] * b[0][0] + a[1][1] * b[1][0]) % n,
         (a[1][0] * b[0][1] + a[1][1] * b[1][1]) % n],
    ]


def _step(k, a, n):
    """The numeric step matrix M(k) = [[2a(2k+1), 0], [2a(2k+1), k+1]]."""
    g = 2 * a * (2 * k + 1) % n
    return [[g, 0], [g, (k + 1) % n]]


def central_sum(t: int, a: int, n: int):
    """``sum_{k<=t} C(2k,k) a^k mod n`` in ``O~(sqrt(t))`` ring operations.

    Tracking numerators ``u_k = k! C(2k,k) a^k`` and ``V_k = k! S_k`` gives the
    division-free recurrence

        [u_{k+1}]   [2a(2k+1)    0 ] [u_k]
        [V_{k+1}] = [2a(2k+1)  k+1 ] [V_k],

    so ``S_t = V_t / t!``.  The product of ``t`` step matrices is split into
    ``c = isqrt(t)`` blocks: one polynomial matrix ``Q(x) = M(x+c-1)...M(x)`` is
    built by a product tree and evaluated at ``x = 0, c, 2c, ...`` -- exactly the
    Bostan-Gaudry-Schost scheme used for ``t! mod n``, one dimension up.

    Returns ``(sum, None)`` or, if ``t!`` shares a factor with ``n``, ``(None, g)``.
    """
    c = isqrt(t)
    if c < 4:
        g = gcd(factorial_mod(t, n), n)
        if g != 1:
            return None, g
        return central_sum_naive(t, a, n), None
    two_a = 2 * a % n
    # leaves M(x + i) as polynomial matrices in x
    leaves = []
    for i in range(c):
        g = [(two_a * (2 * i + 1)) % n, (2 * two_a) % n]   # 2a(2(x+i)+1)
        leaves.append([[g, []], [g, [(i + 1) % n, 1 % n]]])
    # product, later steps on the left: Q = M(x+c-1) ... M(x+1) M(x)
    level = leaves
    while len(level) > 1:
        nxt = []
        for j in range(0, len(level) - 1, 2):
            nxt.append(_mmul(level[j + 1], level[j], n))
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    q_poly = level[0]
    points = [(i * c) % n for i in range(c)]
    evals = [[multipoint_eval(q_poly[r][s] or [0], points, n) for s in range(2)]
             for r in range(2)]
    total = [[1, 0], [0, 1]]
    for i in range(c):
        block = [[evals[r][s][i] for s in range(2)] for r in range(2)]
        total = _nmul(block, total, n)
    for k in range(c * c, t):
        total = _nmul(_step(k, a, n), total, n)
    u_t, v_t = total[0][0] + total[0][1], total[1][0] + total[1][1]
    fact = factorial_mod(t, n)
    g = gcd(fact, n)
    if g != 1:
        return None, g
    return v_t % n * pow(fact, -1, n) % n, None


def central_split(n: int, bases=range(2, 200), stats: dict | None = None):
    """Split ``n`` with the central column: Legendre symbols at ``T = 2^j``.

    For each ``T`` in a geometric sequence, and each base ``a`` with
    ``jacobi(1-4a, n) == -1``, test ``gcd(P_T(a) -+ 1, n)``.  The ``T`` landing in
    ``[(p-1)/2, p-1]`` splits ``n`` through the Legendre symbol; a ``T`` past ``p``
    is caught instead by the factorial's gcd, which is Strassen's mechanism and is
    reported separately in ``stats``.
    """
    if n % 2 == 0:
        return (2, n // 2)
    a_list = [a for a in bases if jacobi(1 - 4 * a, n) == -1][:3]
    t = 2
    while t * t <= 4 * n:
        for a in a_list:
            s, g = central_sum(t, a, n)
            if g is not None:
                if 1 < g < n:
                    if stats is not None:
                        stats.update({"T": t, "mechanism": "factorial (Strassen)"})
                    return (min(g, n // g), max(g, n // g))
                break
            for eps in (1, n - 1):
                g = gcd((s - eps) % n, n)
                if 1 < g < n:
                    if stats is not None:
                        # Classify after the fact.  Inside g's window the value is
                        # the Legendre symbol mod g; below it, s happened to hit
                        # +-1 mod g by chance -- a ~2/g fluke, trial-division-like.
                        in_window = (g - 1) // 2 <= t <= g - 1
                        stats.update({
                            "T": t, "a": a,
                            "mechanism": "Legendre symbol" if in_window
                                         else "coincidence (T below the window)",
                        })
                    return (min(g, n // g), max(g, n // g))
        t *= 2
    return None


def power_series_naive(t: int, a: int, n: int, d: int) -> int:
    """Truncated ``(1 - d^d x)^(-1/d)`` at ``x = a``, mod ``n``, in ``O(t)``.

    The ``d``-th root generalisation of the central column (which is ``d = 2``).
    Coefficient ratio: ``c_{k+1}/c_k = d^(d-1) (1 + d k) / (k + 1)``, and the
    ``d^d`` normalisation keeps every coefficient an integer.

    Modulo a prime ``p`` the truncation at ``T`` in ``[r_d(p), p-1]`` equals
    ``(1 - d^d a)^(r_d(p))`` where ``r_d(p) = (j p - 1)/d`` with
    ``j = p^(-1) mod d``.  When ``d | p-1``, ``j = 1`` and this is the ``d``-th
    power residue character, detected at ``T ~ p/d``.
    """
    s = term = 1 % n
    lead = pow(d, d - 1, n)
    for k in range(t):
        term = term * (lead * (1 + d * k) * a) % n * pow(k + 1, -1, n) % n
        s = (s + term) % n
    return s
