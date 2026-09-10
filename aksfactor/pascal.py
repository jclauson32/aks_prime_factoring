"""Row ``n`` of Pascal's triangle, reduced mod ``n``.

The object this whole package studies is the vector

    R(n) = [ C(n, 0) mod n, C(n, 1) mod n, ..., C(n, n) mod n ]

whose interior is identically zero exactly when ``n`` is prime -- the
truncated ``a = 1`` case of the AKS identity ``(x + a)**n == x**n + a (mod n)``.

Three engines compute entries of ``R(n)``:

``row_entry``    O(1) random access to a single position, for ``n`` of any size.
``row_series``   the whole prefix ``k <= kmax``, via ``(1+x)**n mod (n, x**(kmax+1))``.
``row_exact``    the whole prefix via :func:`math.comb`; a reference, small ``n`` only.
"""

from __future__ import annotations

from math import comb, gcd

from .binom import binom_mod_small

__all__ = [
    "row_entry",
    "row_series",
    "row_exact",
    "row_prefix",
    "row_support",
    "first_nonzero",
    "residue_shape",
]


# --------------------------------------------------------------------------
# Random access
# --------------------------------------------------------------------------

def row_entry(n: int, k: int) -> int:
    """``C(n, k) mod n`` for a single ``k``, without touching any other entry.

    ``n`` may have thousands of digits.  The cost is dominated by the largest
    prime power dividing ``k`` (never by the size of ``n``), because of the
    reduction

        C(n, k) mod n = (n/d) * ( C(n-1, k-1) / (k/d)  mod d ),   d = gcd(n, k)

    proved as Theorem 5 in ``docs/THEORY.md``.  The inner quantity only needs
    ``C(n-1, k-1)`` modulo the *small* number ``k``, which
    :func:`aksfactor.binom.binom_mod_small` supplies via Lucas/Granville.
    """
    if n < 2:
        raise ValueError("row_entry expects n >= 2")
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1 % n
    d = gcd(n, k)
    if d == 1:
        return 0  # Theorem 2: coprime position => n divides C(n, k)
    kp = k // d
    c = binom_mod_small(n - 1, k - 1, k)  # = (k/d) * (m mod d), see docstring
    if c % kp:  # pragma: no cover - guards the Theorem 5 divisibility claim
        raise AssertionError(f"Theorem 5 violated at n={n}, k={k}")
    return (n // d) * (c // kp)


# --------------------------------------------------------------------------
# Whole-prefix engines
# --------------------------------------------------------------------------

def _pack(coeffs: list[int], wbytes: int) -> int:
    return int.from_bytes(
        b"".join(c.to_bytes(wbytes, "little") for c in coeffs), "little"
    )


def _unpack(x: int, wbytes: int, count: int, n: int) -> list[int]:
    raw = x.to_bytes(max((x.bit_length() + 7) // 8, wbytes * count), "little")
    return [
        int.from_bytes(raw[i * wbytes : (i + 1) * wbytes], "little") % n
        for i in range(count)
    ]


def _mul_trunc(a: list[int], b: list[int], kmax: int, n: int) -> list[int]:
    """Truncated product of two polynomials over ``Z/n``, degrees > ``kmax`` dropped."""
    out_len = min(len(a) + len(b) - 1, kmax + 1)
    if min(len(a), len(b)) <= 24:  # schoolbook wins on small operands
        out = [0] * out_len
        for i, ai in enumerate(a):
            if not ai or i > kmax:
                continue
            for j in range(min(len(b), out_len - i)):
                bj = b[j]
                if bj:
                    out[i + j] = (out[i + j] + ai * bj) % n
        return out
    # Kronecker substitution: pack coefficients into one big integer and let
    # CPython's Karatsuba multiplication do the convolution.
    slack = max(len(a), len(b)).bit_length()
    wbytes = (2 * n.bit_length() + slack + 8) // 8 + 1
    prod = _pack(a, wbytes) * _pack(b, wbytes)
    return _unpack(prod, wbytes, out_len, n)


def row_series(n: int, kmax: int) -> list[int]:
    """``[C(n,k) mod n for k in 0..kmax]`` via ``(1+x)**n mod (n, x**(kmax+1))``.

    This is the AKS object itself, truncated: repeated squaring in the ring
    ``(Z/n)[x] / (x**(kmax+1))`` costs ``O(log n)`` polynomial multiplications.
    """
    if n < 2:
        raise ValueError("row_series expects n >= 2")
    kmax = min(kmax, n)
    base = [1 % n, 1 % n][: kmax + 1] or [1 % n]
    result = [1 % n]
    e = n
    while e:
        if e & 1:
            result = _mul_trunc(result, base, kmax, n)
        e >>= 1
        if e:
            base = _mul_trunc(base, base, kmax, n)
    if len(result) < kmax + 1:
        result += [0] * (kmax + 1 - len(result))
    return result[: kmax + 1]


def row_exact(n: int, kmax: int) -> list[int]:
    """``[C(n,k) mod n for k in 0..kmax]`` straight from :func:`math.comb`."""
    kmax = min(kmax, n)
    return [comb(n, k) % n for k in range(kmax + 1)]


def row_prefix(n: int, kmax: int, engine: str = "auto") -> list[int]:
    """Prefix of ``R(n)``.  ``engine`` is ``"auto"``, ``"series"``, ``"exact"``
    or ``"entry"``."""
    if engine == "auto":
        engine = "exact" if n <= 4096 else "series"
    if engine == "exact":
        return row_exact(n, kmax)
    if engine == "series":
        return row_series(n, kmax)
    if engine == "entry":
        return [row_entry(n, k) for k in range(min(kmax, n) + 1)]
    raise ValueError(f"unknown engine {engine!r}")


def row_support(n: int, kmax: int | None = None, engine: str = "auto") -> list[int]:
    """Interior positions ``0 < k < n`` where ``C(n, k) mod n`` is non-zero."""
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    row = row_prefix(n, kmax, engine=engine)
    return [k for k in range(1, kmax + 1) if row[k]]


def first_nonzero(n: int, kmax: int | None = None, engine: str = "auto"):
    """First interior position with a non-zero residue, as ``(k, residue)``.

    Returns ``None`` when the scanned prefix is entirely zero (which, if the
    whole interior was scanned, certifies that ``n`` is prime).
    """
    kmax = n - 1 if kmax is None else min(kmax, n - 1)
    row = row_prefix(n, kmax, engine=engine)
    for k in range(1, kmax + 1):
        if row[k]:
            return k, row[k]
    return None


def residue_shape(n: int, k: int, residue: int | None = None) -> dict:
    """Decompose ``C(n,k) mod n`` as ``x * y`` with ``x = n / gcd(n, k)``.

    This is the decomposition the project set out to test: ``x`` is always a
    proper divisor of ``n`` (Theorem 1), and ``y`` is the leftover cofactor.
    """
    r = row_entry(n, k) if residue is None else residue
    d = gcd(n, k)
    x = n // d
    return {
        "n": n,
        "k": k,
        "residue": r,
        "gcd_nk": d,
        "x": x,
        "y": (r // x) if r and r % x == 0 else None,
        "x_divides_n": n % x == 0,
        "x_divides_residue": (r % x == 0),
        "gcd_residue_n": gcd(r, n),
    }
