"""Round 23: Pascal's triangle over divisibility sequences -- and the one that is ECM.

Replace ``n`` by the ``n``-th term of a *strong divisibility sequence*
(``gcd(a_m, a_n) = a_gcd(m, n)``) and Pascal's triangle generalises:

    [n, k]_a  =  a_n a_(n-1) .. a_(n-k+1) / (a_1 a_2 .. a_k)

is always an integer.  Four instances:

* ``a_n = n``: Pascal's triangle;
* ``a_n = (q^n - 1)/(q - 1)``: the Gaussian (q-) binomials of round 12;
* ``a_n = F_n`` (a Lucas sequence): the Fibonomials;
* ``a_n = W_n``, an elliptic divisibility sequence -- the division polynomials
  of a point ``P`` on an elliptic curve: an *elliptic* Pascal triangle.

Modulo a prime ``p`` each triangle is a Sierpinski gasket, but its base cell is
the *rank of apparition* ``r(p)``, the first ``n`` with ``p | a_n`` (Kummer's
theorem in the mixed radix ``(r, p, p, ...)``, after Knuth and Wilf).  For
Pascal ``r = p``; for the q- and Lucas triangles ``r`` divides ``p - 1`` or
``p +- 1``; for the elliptic triangle ``r`` is the order of ``P`` in
``E(F_p)``, anywhere in the Hasse interval -- and a different curve gives a
different triangle with a different cell.  Jumping to row ``lcm(1..B)`` of
that triangle is exactly stage 1 of Lenstra's ECM.
"""

from __future__ import annotations

from math import gcd


def eds(w1: int, w2: int, w3: int, w4: int, count: int) -> list[int]:
    """``[W_0, .., W_count]`` from the elliptic divisibility recurrence

        W_(2k+1) = W_(k+2) W_k^3 - W_(k-1) W_(k+1)^3,
        W_(2k)   = (W_(k+2) W_(k-1)^2 - W_(k-2) W_(k+1)^2) W_k / W_2.
    """
    w = [0, w1, w2, w3, w4]
    for m in range(5, count + 1):
        k = m // 2
        if m % 2:
            val = w[k + 2] * w[k] ** 3 - w[k - 1] * w[k + 1] ** 3
        else:
            val = (w[k + 2] * w[k - 1] ** 2 - w[k - 2] * w[k + 1] ** 2) * w[k]
            if val % w[2]:
                raise ValueError("not an elliptic divisibility sequence")
            val //= w[2]
        w.append(val)
    return w[: count + 1]


def curve_eds(a: int, b: int, x: int, y: int, count: int) -> list[int]:
    """The EDS of the point ``(x, y)`` on ``y^2 = x^3 + a x + b``: division
    polynomials ``psi_n`` evaluated at the point."""
    if (y * y - x ** 3 - a * x - b) != 0:
        raise ValueError("point not on curve")
    p3 = 3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a * a
    p4 = 4 * y * (x ** 6 + 5 * a * x ** 4 + 20 * b * x ** 3 - 5 * a * a * x * x
                  - 4 * a * b * x - 8 * b * b - a ** 3)
    return eds(1, 2 * y, p3, p4, count)


def fibonacci(count: int) -> list[int]:
    f = [0, 1]
    while len(f) <= count:
        f.append(f[-1] + f[-2])
    return f[: count + 1]


def q_integers(q: int, count: int) -> list[int]:
    return [0] + [(q ** n - 1) // (q - 1) for n in range(1, count + 1)]


def naturals(count: int) -> list[int]:
    return list(range(count + 1))


def nomial(seq, n: int, k: int) -> int:
    """``[n, k]_a`` exactly (raises if the division is not exact)."""
    num = den = 1
    for i in range(1, k + 1):
        num *= seq[n - i + 1]
        den *= seq[i]
    q, r = divmod(num, den)
    if r:
        raise ArithmeticError("generalised binomial is not an integer")
    return q


def rank_of_apparition(seq, p: int) -> int | None:
    """First ``n >= 1`` with ``p | a_n``."""
    for n in range(1, len(seq)):
        if seq[n] % p == 0:
            return n
    return None


def mixed_radix_carries(n: int, k: int, r: int, p: int) -> int:
    """Carries when adding ``k`` and ``n - k`` in the radix ``(r, p, p, ...)``."""
    a, b = k, n - k
    carries = carry = 0
    base = r
    while a or b or carry:
        s = a % base + b % base + carry
        carry = 1 if s >= base else 0
        carries += carry
        a //= base
        b //= base
        base = p
    return carries


def gasket(seq, p: int, rows: int) -> list[str]:
    """Rows ``0..rows-1`` of the triangle mod ``p``: ``#`` non-zero, ``.`` zero."""
    out = []
    for n in range(rows):
        out.append("".join("#" if nomial(seq, n, k) % p else "." for k in range(n + 1)))
    return out


def point_order(a: int, b: int, x: int, y: int, p: int) -> int:
    """Order of ``(x, y)`` on ``y^2 = x^3 + a x + b`` over ``F_p`` (naive)."""
    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return None
        if P == Q:
            lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return x3, (lam * (x1 - x3) - y1) % p

    P0 = (x % p, y % p)
    P, n = P0, 1
    while P is not None:
        P = add(P, P0)
        n += 1
    return n


def eds_mod(w1: int, w2: int, w3: int, w4: int, count: int, m: int) -> list[int]:
    """The EDS modulo ``m`` (needs ``gcd(W_2, m) = 1``)."""
    inv2 = pow(w2 % m, -1, m)
    w = [0, w1 % m, w2 % m, w3 % m, w4 % m]
    for t in range(5, count + 1):
        k = t // 2
        if t % 2:
            val = w[k + 2] * pow(w[k], 3, m) - w[k - 1] * pow(w[k + 1], 3, m)
        else:
            val = (w[k + 2] * w[k - 1] ** 2 - w[k - 2] * w[k + 1] ** 2) * w[k] * inv2
        w.append(val % m)
    return w


def first_shared_row(w1, w2, w3, w4, n: int, limit: int):
    """First ``t`` with ``1 < gcd(W_t, n) < n``: the row where the elliptic
    triangle mod ``n`` first has a whole row divisible by one prime factor."""
    w = eds_mod(w1, w2, w3, w4, limit, n)
    for t in range(1, limit + 1):
        g = gcd(w[t], n)
        if 1 < g < n:
            return t, g
    return None


def eds_at(w1: int, w2: int, w3: int, w4: int, index: int, m: int) -> int:
    """``W_index mod m`` in ``O(log index)`` steps (Shipsey's double-and-add).

    Keeps the block ``W_(c-3) .. W_(c+4)`` around a centre ``c`` and maps it to
    the block around ``2c`` or ``2c + 1`` with the duplication formulas -- the
    elliptic divisibility sequence's own version of the Montgomery ladder, and
    the way to reach row ``lcm(1..B)`` of the elliptic triangle.
    """
    base = eds_mod(w1, w2, w3, w4, 12, m)
    if index <= 12:
        return base[index]
    inv2 = pow(w2 % m, -1, m)
    bits = bin(index)[2:]
    c = int(bits[:3], 2)                      # 4 <= c <= 7
    rest = bits[3:]
    block = {i: base[i] for i in range(c - 3, c + 5)}

    def odd(i):   # W_(2i+1)
        return (block[i + 2] * pow(block[i], 3, m) - block[i - 1] * pow(block[i + 1], 3, m)) % m

    def even(i):  # W_(2i)
        return ((block[i + 2] * block[i - 1] ** 2 - block[i - 2] * block[i + 1] ** 2)
                * block[i] % m * inv2) % m

    for bit in rest:
        new = {}
        lo = 2 * c - 3 + (1 if bit == "1" else 0)
        for t in range(lo, lo + 8):
            i, r = divmod(t, 2)
            new[t] = odd(i) if r else even(i)
        c = 2 * c + (1 if bit == "1" else 0)
        block = new
    return block[c]


def curve_eds_mod(a: int, b: int, x: int, y: int, m: int):
    """``(W_1, W_2, W_3, W_4) mod m`` for the point ``(x, y)`` on ``y^2 = x^3 + ax + b``."""
    p3 = (3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a * a) % m
    p4 = 4 * y * (x ** 6 + 5 * a * x ** 4 + 20 * b * x ** 3 - 5 * a * a * x * x
                  - 4 * a * b * x - 8 * b * b - a ** 3) % m
    return 1, 2 * y % m, p3, p4


def elliptic_triangle_factor(n: int, bound: int, curves: int, rng, stats: dict | None = None):
    """Row ``M = lcm(1..bound)`` of random elliptic triangles mod ``n``.

    Each curve ``y^2 = x^3 + a x + b`` through a random point gives an EDS; its
    triangle mod ``p`` has base cell ``r_p = ord(P mod p)``, and ``W_M`` is
    divisible by ``p`` exactly when ``r_p | M``.  One gcd per curve.
    """
    from .arith import sieve

    m_index = 1
    for ell in sieve(bound):
        e = 1
        while ell ** (e + 1) <= bound:
            e += 1
        m_index *= ell ** e
    for c in range(1, curves + 1):
        a, x, y = rng.randrange(n), rng.randrange(n), rng.randrange(n)
        b = (y * y - x ** 3 - a * x) % n
        w = curve_eds_mod(a, b, x, y, n)
        g = gcd(w[1], n)
        if 1 < g < n:
            return g
        if g == n:
            continue
        g = gcd(eds_at(*w, m_index, n), n)
        if 1 < g < n:
            if stats is not None:
                stats["curves"] = c
            return g
    if stats is not None:
        stats["curves"] = curves
    return None
