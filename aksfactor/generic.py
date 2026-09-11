"""What generic algebra over Z/N achieves, and what every working method adds.

Factoring is the manufacture of a **zero divisor**.  Under CRT,
``Z/N = F_p x F_q``, and a non-trivial ``gcd(x, N)`` means ``x`` vanishes in
exactly one coordinate.  Every method in this repository does the same thing:
build an element of ``Z/N`` and hope it is a zero divisor.

The trouble is that ring operations act **diagonally**: ``+``, ``-`` and ``*``
apply the same rule to both CRT coordinates, so no arrangement of them can tell
``p`` from ``q``.  This module measures the consequence.

``generic_rate``
    the rate at which random straight-line programs over ``Z/N`` -- random base,
    arbitrary ``+ - *``, so deep enough to contain exponentiation by repeated
    squaring -- produce a zero divisor.
``zero_divisor_density``
    ``1 - phi(N)/N``, the rate for an element drawn uniformly at random.

They coincide.  Generic algebra is worth exactly nothing over drawing an element
out of a hat, at any depth.  What every working method adds is a single
**alignment** with a structure attached to ``p``:

* **order alignment** -- choose the exponent to be divisible by ``|G|`` for some
  group ``G`` attached to ``p`` (Pollard ``p-1``, Williams ``p+1``, ECM, class
  groups).  Capped at ``L[1/2]`` by the density of smooth numbers.
* **size alignment** -- choose an interval or window to straddle the magnitude of
  ``p`` (Fermat, Lehman, Strassen, Coppersmith, Harvey).  Capped at ``N**(1/4)``
  by counting, ``N**(1/5)`` with Lehman plus a baby-step/giant-step sweep.

No third alignment is known.  ``docs/FINDINGS.md`` collects the fourteen rounds
that failed to find one.
"""

from __future__ import annotations

import random
from math import gcd

__all__ = [
    "random_program",
    "run_program",
    "zero_divisor_density",
    "generic_rate",
    "alignment_gain",
]


def random_program(depth: int, rng: random.Random) -> list[tuple[int, int, int]]:
    """A random straight-line program: ``depth`` operations drawn from ``+ - *``."""
    return [
        (rng.randrange(3), rng.randrange(0, depth + 3), rng.randrange(0, depth + 3))
        for _ in range(depth)
    ]


def run_program(ops, base: int, n: int) -> int:
    """Evaluate a straight-line program at ``base`` over ``Z/n``.

    Registers start at ``(base, 1, 2)``.  Because ``*`` may square a register,
    depth ``d`` reaches exponents as large as ``2**d`` -- this model *contains*
    Pollard's ``p-1``, which is the point: the difference is not the operations
    available but which exponent gets chosen.
    """
    vals = [base % n, 1 % n, 2 % n]
    for op, i, j in ops:
        x, y = vals[i % len(vals)], vals[j % len(vals)]
        vals.append((x + y) % n if op == 0 else (x - y) % n if op == 1 else x * y % n)
    return vals[-1]


def zero_divisor_density(n: int, p: int, q: int) -> float:
    """``1 - phi(n)/n`` for ``n = pq``: the chance a uniform element splits ``n``."""
    return 1 / p + 1 / q - 1 / n


def generic_rate(n: int, depth: int, trials: int, rng: random.Random) -> float:
    """Measured rate at which random programs of the given depth split ``n``."""
    hits = 0
    for _ in range(trials):
        ops = random_program(depth, rng)
        value = run_program(ops, rng.randrange(2, n), n)
        g = gcd(value, n)
        if 1 < g < n:
            hits += 1
    return hits / trials


def alignment_gain(n: int, bound: int, trials: int, rng: random.Random) -> dict:
    """How much an *aligned* exponent buys over the generic baseline.

    Runs Pollard ``p-1`` -- the same ring operations as a generic program, with
    the single difference that the exponent is ``lcm(1..bound)``, chosen to be
    divisible by the order of the group mod ``p`` whenever that order is smooth.
    """
    from math import lcm

    exponent = 1
    for k in range(2, bound + 1):
        exponent = lcm(exponent, k)
    hits = 0
    for _ in range(trials):
        a = rng.randrange(2, n)
        g = gcd(pow(a, exponent, n) - 1, n)
        if 1 < g < n:
            hits += 1
    return {"bound": bound, "trials": trials, "rate": hits / trials}
