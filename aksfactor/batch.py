"""Round 39: amortisation -- the one exponent that moves.

Every other round in this project asks what a single `N` costs.  This one asks
what a *batch* costs, and the answer is different: the per-number price of
finding smooth parts falls when the numbers are processed together, because a
product tree lets one big multiplication serve the whole batch.

Bernstein's algorithm ("How to find smooth parts of integers"):

1. ``P = prod(primes)``, by a product tree;
2. ``z_i = P mod v_i``, by a remainder tree down the same shape;
3. ``y_i = z_i^(2^e) mod v_i`` with ``2^e >= v_i``, which kills every prime of
   ``v_i`` that divides ``P``, however high its power;
4. ``gcd(y_i, v_i)`` is then the largest divisor of ``v_i`` built from those
   primes -- its smooth part -- and ``v_i`` is smooth exactly when the gcd is
   ``v_i`` itself.

The work is ``O~(M)`` in the total size ``M`` of the batch, against
``O(#primes * #values)`` divisions for trial division.  This is the same
amortisation a sieve performs inside a single factorisation (one pass over an
interval marks every candidate divisible by a prime), which is why the sieving
methods and not the rest are the ones that reach `L[1/2]` and `L[1/3]`.
"""

from __future__ import annotations

from math import gcd


def product_tree(values: list[int]) -> list[list[int]]:
    """Levels of the product tree, leaves first; the last level is one product."""
    if not values:
        return [[1]]
    levels = [list(values)]
    while len(levels[-1]) > 1:
        row = levels[-1]
        levels.append([row[i] * row[i + 1] if i + 1 < len(row) else row[i]
                       for i in range(0, len(row), 2)])
    return levels


def remainder_tree(value: int, levels: list[list[int]]) -> list[int]:
    """``value mod leaf`` for every leaf, walking the tree down once."""
    above = [value % levels[-1][0]]
    for depth in range(len(levels) - 2, -1, -1):
        row = levels[depth]
        below = []
        for i, node in enumerate(row):
            below.append(above[i >> 1] % node)
        above = below
    return above


def smooth_parts(values: list[int], primes: list[int]) -> list[int]:
    """The largest divisor of each value built only from ``primes``."""
    if not values:
        return []
    product = product_tree(primes)[-1][0]
    levels = product_tree(values)
    out = []
    for v, z in zip(values, remainder_tree(product, levels)):
        y = z % v
        for _ in range(v.bit_length().bit_length() + 1):
            y = y * y % v
        out.append(gcd(y, v) if y else v)
    return out


def batch_smooth(values: list[int], primes: list[int]) -> list[bool]:
    """Which values are ``primes``-smooth, all at once."""
    return [part == v for v, part in zip(values, smooth_parts(values, primes))]


def trial_smooth(values: list[int], primes: list[int]) -> list[bool]:
    """The same answer one number at a time -- the thing to beat."""
    out = []
    for v in values:
        m = v
        for ell in primes:
            while m % ell == 0:
                m //= ell
            if m == 1:
                break
        out.append(m == 1)
    return out


def batch_gcd(values: list[int]) -> list[int]:
    """``gcd(v_i, prod_{j != i} v_j)`` for every ``i``, in one pass.

    Bernstein's batch gcd, and the whole of the "Mining your Ps and Qs" attack:
    two RSA moduli that share a prime are both factored by one gcd, and finding
    every such pair in a corpus of ``k`` moduli costs ``O~(k)`` rather than the
    ``k^2 / 2`` gcds the obvious loop takes.

    ``P mod v^2`` is divisible by ``v``, and dividing it out leaves
    ``prod_{j != i} v_j mod v_i``; its gcd with ``v_i`` is what is shared.
    """
    if len(values) < 2:
        return [1] * len(values)
    product = product_tree(values)[-1][0]
    squares = product_tree([v * v for v in values])
    out = []
    for v, rem in zip(values, remainder_tree(product, squares)):
        out.append(gcd(rem // v, v))
    return out
