"""Round 24: Shor's algorithm, simulated classically.

The only known polynomial-time factoring algorithm is quantum.  It uses the
order mechanism -- ``a^r = 1 (mod N)``, then ``gcd(a^(r/2) +- 1, N)`` -- but it
finds ``r`` *directly*, by interference, instead of guessing a multiple of it.
Every classical order method (``p - 1``, ECM, the Pascal and elliptic
triangles) has to blind-exponentiate by ``lcm(1..B)`` and hope ``r`` is smooth.
Shor needs no smoothness at all.

This module runs the algorithm on a classical computer, faithfully: the first
register holds ``Q = 2^t >= N^2`` amplitudes, ``f(x) = a^x mod N`` is evaluated
on all of them, the second register is measured, a quantum Fourier transform
(an FFT here) is applied, the first register is measured, and continued
fractions turn the outcome into ``r``.  The simulation is exact and its cost is
``Theta(Q log Q) = Theta(N^2 log N)`` -- which is the point: the quantum
computer holds the ``Q`` amplitudes for free, and a classical one cannot.
"""

from __future__ import annotations

import cmath
import math
from fractions import Fraction
from math import gcd


def fft(values):
    """Iterative radix-2 FFT (length a power of two)."""
    n = len(values)
    a = list(values)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        w = cmath.exp(-2j * math.pi / length)
        for start in range(0, n, length):
            wk = 1
            half = length // 2
            for k in range(half):
                u = a[start + k]
                v = a[start + k + half] * wk
                a[start + k] = u + v
                a[start + k + half] = u - v
                wk *= w
        length <<= 1
    return a


def shor_distribution(n: int, a: int, rng):
    """Simulate one run up to the final measurement.

    Returns ``(t, probabilities, c)``: the register size, the probability of
    each outcome ``y`` of the first register, and the value ``c = a^x0`` the
    second register collapsed to.
    """
    t = (n * n - 1).bit_length()
    q = 1 << t
    fx = [0] * q
    v = 1
    for x in range(q):
        fx[x] = v
        v = v * a % n
    c = fx[rng.randrange(q)]                  # measuring the second register
    amps = [1.0 if fx[x] == c else 0.0 for x in range(q)]
    norm = math.sqrt(sum(amps))
    spec = fft([z / norm for z in amps])
    probs = [abs(z) ** 2 / q for z in spec]
    return t, probs, c


def order_from_outcome(y: int, t: int, n: int, a: int):
    """Continued fractions of ``y / 2^t``: the smallest denominator ``r < n``
    among the convergents with ``a^r = 1``, trying small multiples too."""
    frac = Fraction(y, 1 << t)
    cf = []
    x = frac
    while True:
        ip = x.numerator // x.denominator
        cf.append(ip)
        rem = x - ip
        if rem == 0:
            break
        x = 1 / rem
    for i in range(1, len(cf) + 1):
        conv = Fraction(cf[i - 1])
        for c in reversed(cf[: i - 1]):
            conv = c + 1 / conv
        d = conv.denominator
        if d >= n:
            break
        for m in range(1, 4):
            if pow(a, d * m, n) == 1:
                return d * m
    return None


def shor(n: int, rng, max_runs: int = 20, stats: dict | None = None,
         lucky_gcd: bool = True):
    """Factor ``n`` (odd, not a prime power) with the simulated algorithm.

    With ``lucky_gcd=False`` a base sharing a factor with ``n`` is redrawn
    instead of returned, so every success goes through the period finding.
    """
    runs = 0
    while runs < max_runs:
        a = rng.randrange(2, n - 1)
        g = gcd(a, n)
        if g > 1:
            if lucky_gcd:
                return g
            continue
        runs += 1
        t, probs, _ = shor_distribution(n, a, rng)
        u, acc, y = rng.random(), 0.0, 0
        for y, pr in enumerate(probs):
            acc += pr
            if acc >= u:
                break
        r = order_from_outcome(y, t, n, a)
        if r is None or r % 2 or pow(a, r // 2, n) == n - 1:
            continue
        for s in (-1, 1):
            g = gcd(pow(a, r // 2, n) + s, n)
            if 1 < g < n:
                if stats is not None:
                    stats.update(runs=runs, t=t, r=r, a=a)
                return g
    if stats is not None:
        stats.update(runs=runs)
    return None


def dft_any(values):
    """DFT of any length by Bluestein's chirp-z transform over ``fft``."""
    n = len(values)
    if n & (n - 1) == 0:
        return fft(values)
    m = 1
    while m < 2 * n - 1:
        m <<= 1
    chirp = [cmath.exp(-1j * math.pi * (k * k % (2 * n)) / n) for k in range(n)]
    a = [values[k] * chirp[k] for k in range(n)] + [0] * (m - n)
    b = [0] * m
    b[0] = chirp[0].conjugate()
    for k in range(1, n):
        b[k] = b[m - k] = chirp[k].conjugate()
    fa, fb = fft(a), fft(b)
    conv = fft([x * y for x, y in zip(fa, fb)])
    conv = [conv[0] / m] + [conv[m - k] / m for k in range(1, m)]   # inverse via reversal
    return [conv[k] * chirp[k] for k in range(n)]


def order_spectrum(n: int, a: int, psi, r: int | None = None):
    """Fourier weights ``|g^(j)|^2`` of ``g(s) = psi(a^s mod n)`` over one period.

    In Shor's algorithm the first register's outcome ``y ~ j Q / r`` carries
    exactly this weight when ``psi`` is the second register's measurement; a
    classical heavy-coefficient search (Goldreich-Levin, Kushilevitz-Mansour)
    for a frequency revealing ``r`` needs about ``1 / max weight`` queries.
    """
    if r is None:
        r, v = 1, a % n
        while v != 1:
            v = v * a % n
            r += 1
    g = [psi(pow(a, s, n)) for s in range(r)]
    spec = dft_any(g)
    return r, [abs(z / r) ** 2 for z in spec]
