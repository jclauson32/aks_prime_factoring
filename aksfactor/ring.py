"""The full AKS object, folded: ``(x + a)**n  mod (n, x**r - 1)``.

Truncating Pascal's row at degree ``kmax`` costs ``O(kmax)`` coefficients.
Folding it instead -- working modulo ``x**r - 1`` -- compresses the *entire*
row into ``r`` numbers for ``O(log n)`` ring multiplications, which is the only
plausible route to reading the row without paying ``Theta(spf(n))``.

It leaks factors.  It just does not leak them fast enough; see
:func:`fold_identity` for exactly why, and ``docs/FINDINGS.md`` for the measured
scaling.
"""

from __future__ import annotations

from math import comb, gcd

from .arith import valuation

__all__ = ["aks_pow", "fold_coefficients", "fold_identity", "fold_attack"]


def _cyclic_mul(a: list[int], b: list[int], r: int, n: int) -> list[int]:
    """Product in ``(Z/n)[x] / (x**r - 1)``."""
    if r <= 24 or min(len(a), len(b)) <= 24:
        out = [0] * r
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        out[(i + j) % r] = (out[(i + j) % r] + ai * bj) % n
        return out
    from .pascal import _pack, _unpack  # local import keeps module graph shallow

    wbytes = (2 * n.bit_length() + r.bit_length() + 8) // 8 + 1
    prod = _pack(a, wbytes) * _pack(b, wbytes)
    flat = _unpack(prod, wbytes, 2 * r - 1, n)
    out = flat[:r]
    for i in range(r, 2 * r - 1):
        out[i - r] = (out[i - r] + flat[i]) % n
    return out


def aks_pow(n: int, r: int, a: int = 1, exponent: int | None = None) -> list[int]:
    """Coefficients of ``(x + a)**n mod (n, x**r - 1)``, index ``j`` first."""
    if r < 1:
        raise ValueError("aks_pow expects r >= 1")
    e = n if exponent is None else exponent
    base = [0] * r
    base[0] = (base[0] + a) % n
    base[1 % r] = (base[1 % r] + 1) % n
    result = [0] * r
    result[0] = 1 % n
    while e:
        if e & 1:
            result = _cyclic_mul(result, base, r, n)
        e >>= 1
        if e:
            base = _cyclic_mul(base, base, r, n)
    return result


def fold_coefficients(n: int, r: int, a: int = 1) -> list[int]:
    """``S_j = sum over k = j (mod r) of C(n,k) * a**(n-k), reduced mod n``."""
    return aks_pow(n, r, a)


def fold_identity(n: int, p: int, r: int, a: int = 1) -> tuple[bool, dict]:
    """Verify why the fold leaks: ``S_j mod p`` is a binomial progression sum.

    Modulo a prime ``p`` with ``p**v || n`` and ``s = n / p**v``, Frobenius gives
    ``(x + a)**n == (x**(p**v) + a)**s (mod p)``, so

        S_j == sum over { i : i * p**v = j (mod r) } of C(s, i) * a**(s-i)   (mod p)

    Each ``S_j`` therefore vanishes mod ``p`` with probability about ``1/p``,
    which is what pins the attack's cost at ``Theta(p)`` coefficients.
    ``s`` must be small enough to enumerate, so this is a small-``n`` check.
    """
    v = valuation(n, p)
    s = n // p**v
    coeffs = fold_coefficients(n, r, a)
    shift = pow(p, v, r)
    for j in range(r):
        want = sum(
            comb(s, i) * pow(a, s - i, p) for i in range(s + 1) if (i * shift) % r == j
        ) % p
        if coeffs[j] % p != want:
            return False, {"n": n, "p": p, "r": r, "a": a, "j": j,
                           "got": coeffs[j] % p, "want": want}
    return True, {"n": n, "p": p, "r": r, "a": a, "s": s}


def fold_attack(
    n: int,
    rmax: int = 200,
    bases: tuple[int, ...] = (1,),
    rmin: int = 2,
) -> dict:
    """Search folds ``(r, a)`` for a coefficient whose gcd with ``n`` is a factor.

    Returns a dict with the factor found (or ``None``) plus ``pairs``: the
    number of ``(r, j, a)`` coefficients inspected.  That count is the honest
    cost measure -- ``docs/FINDINGS.md`` plots it against ``spf(n)``.
    """
    pairs = 0
    for r in range(rmin, rmax + 1):
        for a in bases:
            coeffs = fold_coefficients(n, r, a)
            for j, c in enumerate(coeffs):
                pairs += 1
                g = gcd(c, n)
                if 1 < g < n:
                    return {"factor": g, "cofactor": n // g, "r": r, "j": j,
                            "a": a, "pairs": pairs}
    return {"factor": None, "cofactor": None, "r": None, "j": None,
            "a": None, "pairs": pairs}
