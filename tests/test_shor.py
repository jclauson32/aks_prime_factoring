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
