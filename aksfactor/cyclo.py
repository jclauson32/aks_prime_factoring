"""Escaping the domination argument: factoring in the norm-one subgroup.

``docs/FINDINGS.md`` round 2 argued that generalising Pollard ``p-1`` into the
AKS ring buys nothing, because ``p-1`` divides ``p**k - 1`` for every ``k``, so
``p**k - 1`` is ``B``-smooth only if ``p-1`` already was.  That argument is
correct **about the full unit group**, and it is too strong.

``F_{p^d}*`` has a norm-one subgroup

    ker(N : F_{p^d}* -> F_p*),    of order  (p**d - 1)/(p - 1),

and ``p - 1`` does *not* divide that.  So the smoothness of the relevant group
order genuinely decouples from the smoothness of ``p - 1``, and primes where
``p-1`` is rough while ``(p^d-1)/(p-1)`` is smooth are plentiful.

The catch is landing *inside* the norm-one subgroup without knowing ``p``.  The
trick is to choose the defining polynomial so that its roots have norm one by
construction: a monic ``f`` whose constant term makes the product of the roots
equal ``1``.  For ``d = 2`` that is ``x**2 - a*x + 1``, whose orbit is tracked by
the Lucas sequence ``V_k(a, 1)`` -- i.e. Williams' ``p+1`` method, recovered here
as the ``d = 2`` case.

None of this is polynomial time.  It is a genuine asymmetric statistic, which is
what the first-digit analysis said the framework lacked.
"""

from __future__ import annotations

from math import gcd, isqrt

from .arith import is_prime, sieve

__all__ = [
    "lucas_v",
    "norm_one_search",
    "lucas_v_binomial",
    "williams_pplus1",
    "pollard_pminus1",
    "norm_one_factor",
]


# --------------------------------------------------------------------------
# d = 2: Lucas sequences
# --------------------------------------------------------------------------

def lucas_v(k: int, a: int, n: int) -> int:
    """``V_k(a, 1) mod n``, by the Montgomery ladder.

    ``V`` is the trace of ``z**k`` where ``z`` is a root of ``x**2 - a*x + 1``.
    Because the constant term is ``1``, the two roots multiply to ``1``: ``z``
    has norm one *by construction*, which is the whole point.
    """
    if k == 0:
        return 2 % n
    v0, v1 = 2 % n, a % n
    for bit in bin(k)[2:]:
        if bit == "0":
            v1 = (v0 * v1 - a) % n
            v0 = (v0 * v0 - 2) % n
        else:
            v0 = (v0 * v1 - a) % n
            v1 = (v1 * v1 - 2) % n
    return v0


def lucas_v_binomial(k: int, a: int, n: int) -> int:
    """``V_k(a, 1) mod n`` as an explicit binomial sum.

        V_k = sum_j (-1)**j * (k/(k-j)) * C(k-j, j) * a**(k-2j)

    Slow, and here only to make a point: the statistic that finally breaks the
    symmetry is *still a binomial sum*.  It is simply not a sum along row ``n``.
    """
    from math import comb

    if k == 0:
        return 2 % n
    total = 0
    for j in range(k // 2 + 1):
        term = comb(k - j, j) * k // (k - j)
        total += (-1) ** j * term * pow(a, k - 2 * j)
    return total % n


def _prime_power_ladder(bound: int, exponent_bound: int | None = None):
    """Yield the prime powers whose product is the exponent ``M``.

    ``bound`` caps the *primes*; ``exponent_bound`` separately caps each prime
    *power*.  Keeping them separate matters: ``lcm(1..B)`` only contains ``2**8``
    for ``B = 500``, so a ``p + 1`` of the form ``2**20 * (small primes)`` is
    perfectly ``B``-smooth yet still missed.  Defaulting ``exponent_bound`` well
    above ``bound`` fixes that at negligible cost.
    """
    cap = exponent_bound if exponent_bound is not None else max(bound, 1 << 40)
    for q in sieve(bound):
        power = q
        while power * q <= cap:
            power *= q
        yield power


def williams_pplus1(n: int, bound: int = 10000, base: int = 5,
                    exponent_bound: int | None = None) -> int | None:
    """Williams ``p+1``: the ``d = 2`` norm-one method.

    Finds a prime factor ``p`` when ``p + 1`` is ``bound``-smooth *and* the
    discriminant ``base**2 - 4`` is a non-residue mod ``p`` (otherwise the root
    falls into ``F_p`` and the method degenerates to ``p-1``).  Returns a
    non-trivial divisor or ``None``.
    """
    v = base % n
    for power in _prime_power_ladder(bound, exponent_bound):
        v = lucas_v(power, v, n)
        g = gcd(v - 2, n)
        if 1 < g < n:
            return g
        if g == n:
            return None
    g = gcd(v - 2, n)
    return g if 1 < g < n else None


def pollard_pminus1(n: int, bound: int = 10000, base: int = 2,
                    exponent_bound: int | None = None) -> int | None:
    """Classical Pollard ``p-1``, for side-by-side comparison."""
    a = base % n
    for power in _prime_power_ladder(bound, exponent_bound):
        a = pow(a, power, n)
        g = gcd(a - 1, n)
        if 1 < g < n:
            return g
        if g == n:
            return None
    g = gcd(a - 1, n)
    return g if 1 < g < n else None


# --------------------------------------------------------------------------
# general d: norm-one elements of a degree-d quotient
# --------------------------------------------------------------------------

def _polymod(u: list[int], f: list[int], n: int) -> list[int]:
    """Reduce ``u`` modulo the monic polynomial ``f``, over ``Z/n``."""
    d = len(f) - 1
    u = u[:]
    for i in range(len(u) - 1, d - 1, -1):
        c = u[i]
        if c:
            u[i] = 0
            for j in range(d):
                u[i - d + j] = (u[i - d + j] - c * f[j]) % n
    out = u[:d]
    return out + [0] * (d - len(out))


def _mulmod(u: list[int], v: list[int], f: list[int], n: int) -> list[int]:
    prod = [0] * (len(u) + len(v) - 1)
    for i, ui in enumerate(u):
        if ui:
            for j, vj in enumerate(v):
                if vj:
                    prod[i + j] = (prod[i + j] + ui * vj) % n
    return _polymod(prod, f, n)


def _powmod(u: list[int], e: int, f: list[int], n: int) -> list[int]:
    d = len(f) - 1
    result = [0] * d
    result[0] = 1 % n
    base = u[:]
    while e:
        if e & 1:
            result = _mulmod(result, base, f, n)
        e >>= 1
        if e:
            base = _mulmod(base, base, f, n)
    return result


def _det(matrix: list[list[int]]) -> int:
    """Fraction-free (Bareiss) determinant over the integers."""
    m = [row[:] for row in matrix]
    size = len(m)
    prev = 1
    sign = 1
    for k in range(size - 1):
        if m[k][k] == 0:
            for s in range(k + 1, size):
                if m[s][k] != 0:
                    m[k], m[s] = m[s], m[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
        prev = m[k][k]
    return sign * m[size - 1][size - 1]


def _norm(w: list[int], f: list[int], n: int) -> int:
    """Norm of ``w`` from ``(Z/n)[x]/(f)`` down to ``Z/n``.

    Equals the resultant ``Res(f, w)``, so it vanishes mod ``p`` exactly when
    ``w`` is a zero divisor in some component of the reduction mod ``p`` -- i.e.
    when *any* component satisfies the relation, not all of them.
    """
    d = len(f) - 1
    cols = []
    cur = w[:] + [0] * (d - len(w))
    for _ in range(d):
        cols.append(cur[:])
        cur = _mulmod(cur, [0, 1 % n], f, n)
    matrix = [[cols[j][i] for j in range(d)] for i in range(d)]
    return _det(matrix) % n


def norm_one_factor(
    n: int,
    degree: int = 2,
    bound: int = 10000,
    coeffs: tuple[int, ...] | None = None,
    exponent_bound: int | None = None,
) -> int | None:
    """Norm-one cyclotomic method in a degree-``degree`` quotient.

    Builds ``f = x**d - c_{d-1} x**(d-1) - ... - c_1 x + (-1)**d`` so the roots
    multiply to ``1``, raises ``x`` to ``lcm(1..bound)`` stagewise, and takes
    ``gcd(Norm(z**M - 1), n)``.  Succeeds on a prime ``p`` when the order of the
    root divides ``M`` in some component -- a condition on ``(p**d - 1)/(p - 1)``,
    not on ``p - 1``.
    """
    if degree < 2:
        raise ValueError("norm_one_factor expects degree >= 2")
    if coeffs is None:
        coeffs = tuple(range(1, degree))
    if len(coeffs) != degree - 1:
        raise ValueError("coeffs must have length degree-1")
    f = [(-1) ** degree % n] + [(-c) % n for c in coeffs] + [1 % n]
    z = [0] * degree
    z[1] = 1 % n
    for power in _prime_power_ladder(bound, exponent_bound):
        z = _powmod(z, power, f, n)
        w = z[:]
        w[0] = (w[0] - 1) % n
        g = gcd(_norm(w, f, n), n)
        if 1 < g < n:
            return g
        if g == n:
            return None
    return None


def norm_one_search(
    n: int,
    bound: int = 10000,
    degree: int = 2,
    bases=range(3, 40),
    exponent_bound: int | None = None,
) -> tuple[int, int] | None:
    """Try several bases until the norm-one method splits ``n``.

    Retries are not optional.  For ``degree == 2`` the root only lands in
    ``F_{p^2}`` when the discriminant ``base**2 - 4`` is a non-residue mod ``p``;
    otherwise it falls back into ``F_p`` and the method silently degenerates to
    ``p-1``.  That is a coin flip per base, so a single base failing says nothing.

    Returns ``(factor, base)`` or ``None``.
    """
    for base in bases:
        if degree == 2:
            got = williams_pplus1(n, bound=bound, base=base,
                                  exponent_bound=exponent_bound)
        else:
            coeffs = tuple((base + i) for i in range(degree - 1))
            got = norm_one_factor(n, degree=degree, bound=bound, coeffs=coeffs,
                                  exponent_bound=exponent_bound)
        if got:
            return got, base
    return None
