import random

from aksfactor.shor import fft, order_from_outcome, shor, shor_distribution


def test_fft_matches_dft():
    import cmath

    rng = random.Random(1)
    xs = [complex(rng.random(), rng.random()) for _ in range(16)]
    got = fft(xs)
    for k in range(16):
        want = sum(x * cmath.exp(-2j * cmath.pi * j * k / 16) for j, x in enumerate(xs))
        assert abs(got[k] - want) < 1e-9


def test_peaks_sit_at_multiples_of_q_over_r():
    t, probs, _ = shor_distribution(21, 2, random.Random(3))
    q = 1 << t
    top = sorted(range(q), key=lambda y: -probs[y])[:6]
    assert all(abs(y * 6 / q - round(y * 6 / q)) < 0.01 for y in top)
    assert any(order_from_outcome(y, t, 21, 2) == 6 for y in top)


def test_shor_factors_through_period_finding():
    rng = random.Random(2)
    for n in (15, 21, 35, 77, 143):
        st = {}
        g = shor(n, rng, stats=st, lucky_gcd=False)
        assert g and 1 < g < n and n % g == 0
        assert pow(st["a"], st["r"], n) == 1


def test_bluestein_and_flat_spectra():
    import cmath
    import math
    from math import gcd

    from aksfactor.shor import dft_any, order_spectrum

    rng = random.Random(4)
    xs = [complex(rng.random(), rng.random()) for _ in range(45)]
    got = dft_any(xs)
    want = [sum(x * cmath.exp(-2j * math.pi * j * k / 45) for j, x in enumerate(xs)) for k in range(45)]
    assert max(abs(g - w) for g, w in zip(got, want)) < 1e-9
    n, a = 101 * 103, 5
    r, w = order_spectrum(n, a, lambda v: cmath.exp(2j * math.pi * 17 * v / n))
    assert abs(sum(w) - 1) < 1e-9
    assert max(w[j] for j in range(1, r) if gcd(j, r) == 1) * r < 30
