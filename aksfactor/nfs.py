"""Round 26: a toy number field sieve -- the sign mechanism at L[1/3].

The quadratic sieve asks numbers of size ~ sqrt(N) to be smooth.  The number
field sieve asks two numbers of size ~ N^(1/(d+1)) each to be smooth, by
working in two rings at once: ``Z`` and ``Z[alpha]`` with ``f(alpha) = 0`` and
``f(m) = 0 (mod N)``.  The ring map ``phi: alpha -> m`` sends a square in
``Z[alpha]`` to a square mod ``N`` -- the congruence of squares again, with the
square root taken in a number field.

Pieces (all pure Python, all small):

* base-``m`` polynomial selection, ``f`` monic of degree ``d``;
* a line sieve over ``(a, b)`` on ``a + b m`` and the norm ``F(a, -b)``;
* rational and algebraic factor bases, quadratic characters for the
  obstructions a GF(2) dependency does not see;
* GF(2) elimination;
* the algebraic square root by Newton lifting from an inert prime.
"""

from __future__ import annotations

from math import gcd, isqrt, log

from .arith import is_prime, sieve


# ---------------------------------------------------------------- polynomials
def poly_mulmod(a, b, f, mod=None):
    """``a * b`` in ``Z[x]/(f)`` (``f`` monic, coefficient lists low-first)."""
    d = len(f) - 1
    prod = [0] * (2 * d - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                prod[i + j] += x * y
    for k in range(len(prod) - 1, d - 1, -1):
        c = prod[k]
        if c:
            for i in range(d + 1):
                prod[k - d + i] -= c * f[i]
    out = prod[:d]
    if mod is not None:
        out = [x % mod for x in out]
    return out


def poly_eval(coeffs, x, mod=None):
    v = 0
    for c in reversed(coeffs):
        v = v * x + c
        if mod is not None:
            v %= mod
    return v


def poly_powmod(a, e, f, mod):
    result = [1] + [0] * (len(f) - 2)
    base = [x % mod for x in a]
    while e:
        if e & 1:
            result = poly_mulmod(result, base, f, mod)
        base = poly_mulmod(base, base, f, mod)
        e >>= 1
    return result


def roots_mod(f, p):
    """Roots of ``f`` modulo the prime ``p`` (brute force; ``p`` is small)."""
    return [r for r in range(p) if poly_eval(f, r, p) == 0]


def _base_m(n: int, m: int, d: int):
    coeffs, rest = [], n
    for _ in range(d):
        coeffs.append(rest % m)
        rest //= m
    coeffs.append(rest)
    return coeffs


def _has_rational_root(coeffs) -> bool:
    c0 = abs(coeffs[0])
    if c0 == 0:
        return True
    for t in range(1, isqrt(c0) + 1):
        if c0 % t == 0:
            for r in (t, c0 // t):
                if poly_eval(coeffs, r) == 0 or poly_eval(coeffs, -r) == 0:
                    return True
    return False


def root_score(f, bound: int = 100) -> float:
    """A crude Murphy-alpha: ``sum (roots mod p) log p / (p - 1)`` over small
    primes.  More roots mod small primes means norms divisible by small primes
    more often -- smaller effective norms."""
    return sum(len(roots_mod(f, p)) * log(p) / (p - 1) for p in sieve(bound))


def select_polynomial(n: int, d: int = 3, tries: int = 1):
    """Base-``m`` expansions ``f(m) = n`` with ``f`` monic and irreducible.

    With ``tries = 1`` returns ``m = floor(n^(1/d))``; otherwise ranks up to
    ``tries`` candidates ``m`` just below it by ``root_score`` and returns the
    list best-first as ``[(m, f), ...]``.
    """
    m0 = int(round(n ** (1.0 / d)))
    while m0 ** d > n:
        m0 -= 1
    while (m0 + 1) ** d <= n:
        m0 += 1
    found = []
    m = m0
    while len(found) < tries and m > 1:
        f = _base_m(n, m, d)
        if f[-1] == 1 and not _has_rational_root(f):
            found.append((m, f))
        m -= 1
        if f[-1] != 1:
            break
    if not found:
        raise ValueError("no monic irreducible base-m polynomial")
    if tries == 1:
        return found[0]
    return sorted(found, key=lambda mf: -root_score(mf[1]))


def norm(f, a: int, b: int) -> int:
    """``N(a - b alpha) = b^d f(a/b)`` for monic ``f``, homogenised."""
    d = len(f) - 1
    return sum(c * a ** i * b ** (d - i) for i, c in enumerate(f))


# ---------------------------------------------------------------- the sieve
def factor_bases(f, rational_bound: int, algebraic_bound: int, characters: int):
    rat = sieve(rational_bound)
    alg = [(p, r) for p in sieve(algebraic_bound) for r in roots_mod(f, p)]
    df = [i * c for i, c in enumerate(f)][1:]
    chars = []
    q = algebraic_bound + 1
    while len(chars) < characters:
        q += 1
        if not is_prime(q):
            continue
        for s in roots_mod(f, q):
            if poly_eval(df, s, q) % q:
                chars.append((q, s))
                if len(chars) == characters:
                    break
    return rat, alg, chars


def _trial(v, primes):
    e = []
    for p in primes:
        k = 0
        while v % p == 0:
            v //= p
            k += 1
        e.append(k)
    return e, v


def collect_relations(n, m, f, rat, alg, need, a_width=4000, b_max=4000, stats=None):
    """Line sieve: for each ``b``, ``a`` in ``[-a_width, a_width]``, keep
    ``(a, b)`` with ``a + b m`` smooth over ``rat`` and ``N(a + b alpha)``... over
    ``alg``.  We use the element ``a - b alpha`` and rational side ``a - b m``."""
    rels = []
    alg_primes = sorted({p for p, _ in alg})
    roots_by_p = {}
    for p, r in alg:
        roots_by_p.setdefault(p, []).append(r)
    rlog = {p: log(p) for p in rat}
    alog = {p: log(p) for p in alg_primes}
    size = 2 * a_width + 1
    tried = 0
    for b in range(1, b_max + 1):
        # rational side: a - b m = 0 mod p  <=>  a = b m mod p
        racc = [0.0] * size
        for p in rat:
            start = (b * m + a_width) % p
            lg = rlog[p]
            for i in range(start, size, p):
                racc[i] += lg
        # algebraic side: a - b alpha has norm divisible by p when a = b r mod p
        aacc = [0.0] * size
        for p in alg_primes:
            lg = alog[p]
            for r in roots_by_p[p]:
                start = (b * r + a_width) % p
                for i in range(start, size, p):
                    aacc[i] += lg
        rthr = log(b * m + a_width) - 2.5 * log(rat[-1])
        athr = 2.5 * log(alg_primes[-1])
        for i in range(size):
            if racc[i] < rthr:
                continue
            a = i - a_width
            if a == 0 or gcd(a, b) != 1:
                continue
            nv = norm(f, a, b)
            if nv == 0 or aacc[i] < log(abs(nv)) - athr:
                continue
            rv = a - b * m
            tried += 1
            re, rrest = _trial(abs(rv), rat)
            if rrest != 1:
                continue
            _, arest = _trial(abs(nv), alg_primes)
            if arest != 1:
                continue
            rels.append((a, b, re, rv < 0))
            if len(rels) >= need:
                if stats is not None:
                    stats.update(b=b, candidates=tried)
                return rels
    if stats is not None:
        stats.update(b=b_max, candidates=tried)
    return rels


def relation_vector(a, b, re, rneg, f, rat, alg, chars):
    """GF(2) vector: sign and rational exponents, algebraic ideal exponents
    ``(p, r)`` with ``r = a/b mod p``, quadratic characters of ``a - b s mod q``."""
    bits = [1 if rneg else 0] + [e & 1 for e in re]
    nv = abs(norm(f, a, b))
    for p, r in alg:
        if (a - b * r) % p:
            bits.append(0)
            continue
        k = 0
        while nv % p == 0:
            nv //= p
            k += 1
        bits.append(k & 1)
    for q, s in chars:
        bits.append(0 if pow((a - b * s) % q, (q - 1) // 2, q) == 1 else 1)
    bits.append(1 if norm(f, a, b) < 0 else 0)
    return bits


def dependencies(vectors):
    """All GF(2) dependencies among ``vectors``, as bitmasks over rows."""
    width = len(vectors[0])
    pivots = {}
    deps = []
    for i, v in enumerate(vectors):
        mask = sum(1 << j for j, x in enumerate(v) if x)
        comb = 1 << i
        for bit in range(width):
            if not (mask >> bit) & 1:
                continue
            if bit in pivots:
                pm, pc = pivots[bit]
                mask ^= pm
                comb ^= pc
            else:
                pivots[bit] = (mask, comb)
                break
        if mask == 0:
            deps.append(comb)
    return deps


# ---------------------------------------------------------------- square root
def inert_prime(f, start: int = 1000):
    """A prime ``p`` with ``f`` irreducible mod ``p`` (for degree 3: no roots)."""
    p = start
    while True:
        p += 1
        if is_prime(p) and not roots_mod(f, p) and len(f) - 1 == 3:
            return p


def _field_sqrt(g, f, p):
    """A square root of ``g`` in ``F_p[x]/(f)``, ``f`` irreducible of degree 3.

    ``q = p^3``; use Tonelli-Shanks in the cyclic group ``F_q*``."""
    q = p ** 3
    if poly_powmod(g, (q - 1) // 2, f, p) != [1, 0, 0]:
        return None
    s, t = 0, q - 1
    while t % 2 == 0:
        t //= 2
        s += 1
    # a non-residue
    z = [2, 1, 0]
    while poly_powmod(z, (q - 1) // 2, f, p) == [1, 0, 0]:
        z = [z[0] + 1, 1, 0]
    c = poly_powmod(z, t, f, p)
    x = poly_powmod(g, (t + 1) // 2, f, p)
    b = poly_powmod(g, t, f, p)
    mm = s
    one = [1, 0, 0]
    while b != one:
        i, bb = 0, b
        while bb != one:
            bb = poly_mulmod(bb, bb, f, p)
            i += 1
        cc = c
        for _ in range(mm - i - 1):
            cc = poly_mulmod(cc, cc, f, p)
        x = poly_mulmod(x, cc, f, p)
        c = poly_mulmod(cc, cc, f, p)
        b = poly_mulmod(b, c, f, p)
        mm = i
    return x


def _inverse_mod_pk(g, f, p, k):
    """Inverse of ``g`` in ``(Z/p^k)[x]/(f)`` by Newton from ``F_p``."""
    q = p ** 3
    inv = poly_powmod(g, q - 2, f, p)
    mod = p
    while mod < p ** k:
        mod = min(mod * mod, p ** k)
        e = poly_mulmod(g, inv, f, mod)
        two_minus = [(-x) % mod for x in e]
        two_minus[0] = (two_minus[0] + 2) % mod
        inv = poly_mulmod(inv, two_minus, f, mod)
    return inv


def algebraic_sqrt(gamma, f, bound_bits: int):
    """``beta`` with ``beta^2 = gamma`` in ``Z[alpha]``, coefficients below
    ``2^bound_bits``: square root mod an inert prime, lifted by Newton's
    iteration ``beta <- (beta + gamma / beta) / 2``, doubling the ``p``-adic
    precision each step until it exceeds ``2^(bound_bits+2)``."""
    p = inert_prime(f)
    beta = _field_sqrt([x % p for x in gamma], f, p)
    if beta is None:
        return None
    target = 1 << (bound_bits + 2)
    e, mod = 1, p
    while mod < target:
        e, mod = 2 * e, mod * mod                 # track the exponent exactly
        inv = _inverse_mod_pk(beta, f, p, e)
        t = poly_mulmod([x % mod for x in gamma], inv, f, mod)
        half = pow(2, -1, mod)
        beta = [((x + y) * half) % mod for x, y in zip(beta, t)]
    half_mod = mod // 2
    beta = [x - mod if x > half_mod else x for x in beta]
    if poly_mulmod(beta, beta, f) != list(gamma):
        return None
    return beta


# ---------------------------------------------------------------- the whole thing
def number_field_sieve(n: int, d: int = 3, rational_bound: int = 2000,
                       algebraic_bound: int = 2000, characters: int = 20,
                       a_width: int = 4000, b_max: int = 4000, stats: dict | None = None,
                       polynomials: int = 3):
    """A factor of ``n`` (odd composite, no small factors) or ``None``.

    Tries the base-``m`` polynomials best-first by ``root_score`` until one
    yields enough relations and a split."""
    for m, f in select_polynomial(n, d, tries=12)[:polynomials]:
        g = _nfs_with(n, d, m, f, rational_bound, algebraic_bound, characters,
                      a_width, b_max, stats)
        if g:
            return g
    return None


def _nfs_with(n, d, m, f, rational_bound, algebraic_bound, characters, a_width, b_max, stats):
    rat, alg, chars = factor_bases(f, rational_bound, algebraic_bound, characters)
    width = 2 + len(rat) + len(alg) + len(chars)
    need = width + 10
    st = {}
    rels = collect_relations(n, m, f, rat, alg, need, a_width, b_max, st)
    if stats is not None:
        stats.update(m=m, f=f, width=width, relations=len(rels), **st)
    if len(rels) <= width:
        return None
    vecs = [relation_vector(a, b, re, rneg, f, rat, alg, chars) for a, b, re, rneg in rels]
    df = [i * c for i, c in enumerate(f)][1:] + [0]
    df = df[: len(f) - 1]
    fprime_m = poly_eval(df, m, n)
    for comb in dependencies(vecs):
        chosen = [rels[i] for i in range(len(rels)) if (comb >> i) & 1]
        # rational side: product of (a - b m) is a square of an integer
        total = [0] * len(rat)
        neg = 0
        for a, b, re, rneg in chosen:
            total = [x + y for x, y in zip(total, re)]
            neg ^= rneg
        if neg or any(x % 2 for x in total):
            continue
        x_rat = 1
        for p, e in zip(rat, total):
            x_rat = x_rat * pow(p, e // 2, n) % n
        # algebraic side: gamma = f'(alpha)^2 * prod (a - b alpha)
        gamma = [1] + [0] * (d - 1)
        for a, b, _, _ in chosen:
            gamma = poly_mulmod(gamma, [a, -b] + [0] * (d - 2), f)
        gamma = poly_mulmod(gamma, poly_mulmod(df, df, f), f)
        bound = max(abs(c) for c in gamma).bit_length() + 16
        beta = algebraic_sqrt(gamma, f, bound)
        if beta is None:
            continue
        y_alg = poly_eval(beta, m, n)
        for s in (1, -1):
            g = gcd((y_alg - s * x_rat * fprime_m) % n, n)
            if 1 < g < n:
                return g
    return None
