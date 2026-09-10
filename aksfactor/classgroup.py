"""Class groups of imaginary quadratic orders: a group order that is *dense*.

Round 4 ended on a sharpened requirement.  Every construction in this repository
succeeds when some integer attached to ``p`` is smooth, and the integers on offer
have all come from sparse, structured sets:

* the classical Pascal row: the period is ``p`` -- one value, no choice;
* the q-row: the period is ``ord_p(q)`` -- the divisor lattice of ``p - 1``;
* the norm-one subgroup: the order is ``Phi_d(p)`` -- one value per ``(p, d)``.

What ECM has instead is an *interval*: ``#E(F_p) = p + 1 - t`` with
``|t| <= 2 sqrt(p)``, so a fresh curve is a genuinely fresh number and one can
simply keep sampling until a smooth order appears.

Class groups give the same thing without elliptic curves.  For a discriminant
``D = -kn`` the class number ``h(D)`` sits near ``sqrt(|D|)`` and moves
essentially arbitrarily with ``k`` -- a dense family, indexed by a parameter we
choose, computable without knowing ``p``.  And ambiguous forms (elements of order
dividing 2) *are* factorizations of ``D``, so a smooth ``h(D)`` hands back a
factor of ``n``.

This is the Schnorr-Lenstra method.  It is implemented here because it is the
one construction that meets round 4's requirement, and having it in the same
codebase makes the contrast measurable rather than asserted.

Forms are triples ``(a, b, c)`` of discriminant ``D = b*b - 4*a*c < 0``, always
kept reduced.
"""

from __future__ import annotations

from math import gcd, isqrt

from .arith import sieve

__all__ = [
    "reduce_form",
    "identity_form",
    "compose",
    "square_form",
    "pow_form",
    "class_number_bruteforce",
    "ambiguous_factor",
    "schnorr_lenstra_split",
]


def _xgcd(a: int, b: int) -> tuple[int, int, int]:
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def reduce_form(a: int, b: int, c: int) -> tuple[int, int, int]:
    """The unique reduced form equivalent to ``(a, b, c)``.

    Reduced means ``-a < b <= a <= c``, and ``b >= 0`` when ``a == c``.
    """
    if a <= 0:
        raise ValueError("reduce_form expects a positive definite form")
    d = b * b - 4 * a * c
    while True:
        # normalise b into (-a, a]
        if not (-a < b <= a):
            q = (a - b) // (2 * a)
            b += 2 * a * q
            c = (b * b - d) // (4 * a)
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        break
    return a, b, c


def identity_form(d: int) -> tuple[int, int, int]:
    """The principal form of discriminant ``d``."""
    b = d % 2
    return reduce_form(1, b, (b * b - d) // 4)


def _represents_coprime(f: tuple[int, int, int], m: int, limit: int = 200):
    """Find ``(x, y)`` coprime with ``gcd(f(x,y), m) == 1``, plus a unimodular completion."""
    a, b, c = f
    for bound in range(1, limit):
        for x in range(-bound, bound + 1):
            for y in range(-bound, bound + 1):
                if gcd(x, y) != 1:
                    continue
                val = a * x * x + b * x * y + c * y * y
                if val > 0 and gcd(val, m) == 1:
                    g, w, z = _xgcd(x, y)  # w*x + z*y = g
                    if g < 0:
                        # _xgcd can return a negative gcd for negative inputs;
                        # left unnormalised this yields det = -1, an *improper*
                        # equivalence that lands in a different class.
                        g, w, z = -g, -w, -z
                    assert g == 1 and x * w + y * z == 1
                    assert x * w - y * (-z) == 1  # det of [[x, -z], [y, w]]
                    return x, y, -z, w
    raise RuntimeError("no coprime representation found")


def compose(f1: tuple[int, int, int], f2: tuple[int, int, int], d: int):
    """Gauss composition of two forms of discriminant ``d``.

    Works by first replacing ``f1`` with an equivalent form whose leading
    coefficient is coprime to ``a2``.  In that case the composite is simply

        A = a1' * a2,   B = CRT(b1' mod 2a1', b2 mod 2a2),   C = (B*B - d)/(4A)

    which needs no case analysis -- the awkward part of the textbook algorithm is
    exactly the non-coprime case, and this sidesteps it.
    """
    a1, b1, c1 = reduce_form(*f1)
    a2, b2, c2 = reduce_form(*f2)
    if gcd(a1, a2) != 1:
        x, y, z, w = _represents_coprime((a1, b1, c1), a2)
        na = a1 * x * x + b1 * x * y + c1 * y * y
        nb = 2 * (a1 * x * z + c1 * y * w) + b1 * (x * w + y * z)
        a1, b1 = na, nb
        c1 = (b1 * b1 - d) // (4 * a1)
    # CRT: B = b1 + 2*a1*t  with  a1*t == (b2-b1)/2  (mod a2)
    t = ((b2 - b1) // 2) * pow(a1 % a2, -1, a2) % a2 if a2 != 1 else 0
    b = b1 + 2 * a1 * t
    a = a1 * a2
    c = (b * b - d) // (4 * a)
    return reduce_form(a, b, c)


def square_form(f: tuple[int, int, int], d: int):
    """``f`` composed with itself."""
    return compose(f, f, d)


def pow_form(f: tuple[int, int, int], e: int, d: int):
    """``f**e`` by square-and-multiply."""
    result = identity_form(d)
    base = reduce_form(*f)
    while e:
        if e & 1:
            result = compose(result, base, d)
        e >>= 1
        if e:
            base = square_form(base, d)
    return result


def class_number_bruteforce(d: int) -> int:
    """``h(d)`` by enumerating reduced forms.  Only for small ``|d|``."""
    if d >= 0 or d % 4 not in (0, 1):
        raise ValueError("discriminant must be negative and 0 or 1 mod 4")
    count = 0
    bound = isqrt(-d // 3) + 1
    for a in range(1, bound + 1):
        for b in range(-a + 1, a + 1):
            num = b * b - d
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if gcd(gcd(a, b), c) != 1:
                continue  # only primitive forms
            if a == c and b < 0:
                continue
            count += 1
    return count


def ambiguous_factor(f: tuple[int, int, int], d: int, n: int) -> int | None:
    """Turn an ambiguous form into a factor of ``n``, if it yields one.

    A reduced form of order dividing 2 has ``b == 0``, ``a == b`` or ``a == c``,
    and each case exhibits a factorization of the discriminant.
    """
    a, b, c = f
    candidates = []
    if b == 0:
        candidates += [a, c]
    if a == b:
        candidates += [a, a - 4 * c]
    if a == c:
        candidates += [b - 2 * a, b + 2 * a]
    for value in candidates:
        g = gcd(abs(value), n)
        if 1 < g < n:
            return g
    return None


def schnorr_lenstra_split(
    n: int,
    multipliers=range(1, 40),
    bound: int = 500,
    forms_per_disc: int = 3,
) -> tuple[int, int] | None:
    """Split ``n`` via ambiguous forms in ``Cl(-k*n)``, sweeping the multiplier ``k``.

    For each ``k`` the group order ``h(-kn)`` is a fresh number near
    ``sqrt(kn)``; when its odd part is ``bound``-smooth, raising a random form to
    the smooth exponent leaves an element of 2-power order, and squaring down to
    the identity passes through an ambiguous form -- which factors ``-kn``, and
    with luck ``n``.

    Returns ``(factor, multiplier)`` or ``None``.
    """
    odd_ladder = []
    for q in sieve(bound):
        if q == 2:
            continue
        power = q
        while power * q <= bound:
            power *= q
        odd_ladder.append(power)

    for k in multipliers:
        d = -k * n
        if d % 4 not in (0, 1):
            continue
        for seed in range(1, forms_per_disc + 1):
            a = seed + 1
            # find a form (a, b, *) of discriminant d
            b = None
            for cand in range(d % 2, 2 * a, 2):
                if (cand * cand - d) % (4 * a) == 0:
                    b = cand
                    break
            if b is None:
                continue
            f = reduce_form(a, b, (b * b - d) // (4 * a))
            try:
                for power in odd_ladder:
                    f = pow_form(f, power, d)
                # now squeeze out the 2-part
                identity = identity_form(d)
                for _ in range(64):
                    if f == identity:
                        break
                    nxt = square_form(f, d)
                    if nxt == identity:
                        got = ambiguous_factor(f, d, n)
                        if got:
                            return got, k
                        break
                    f = nxt
            except (RuntimeError, ValueError, ZeroDivisionError):
                continue
    return None
