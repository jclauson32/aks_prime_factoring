import random

from aksfactor.evolve import (
    as_polynomial,
    cost,
    evolve,
    exponent_of,
    optimised_pminus1,
    pminus1_success,
    run,
    score,
)


def test_program_semantics():
    prog = [("pow", 0, 3), ("mul", 2, 0), ("sub", 3, 1)]      # a^3 * a - 1
    assert run(prog, 1000003, 5) == (5 ** 4 - 1) % 1000003
    assert as_polynomial(prog) == {4: 1, 0: -1}
    assert exponent_of(prog[:2]) == 4
    assert cost(prog) == (1.5849625007211563 + 1) + 1


def test_pminus1_success_counts_proper_splits():
    n = 1009 * 1013                          # 1008 = 2^4 * 3^2 * 7, 1012 = 2^2 * 11 * 23
    cases = [(n, a) for a in range(2, 40)]
    assert pminus1_success(2 ** 4 * 3 ** 2 * 7, cases) > 0.9
    assert pminus1_success(1, cases) == 0.0


def test_selection_beats_random():
    rng = random.Random(1)
    primes = [p for p in range(1025, 2048) if all(p % d for d in range(2, 46))]
    cases = []
    while len(cases) < 150:
        p, q = rng.sample(primes, 2)
        cases.append((p * q, rng.randrange(2, p * q - 1)))
    best = evolve(cases, 12, 30, 25, 40, rng)
    e, _ = optimised_pminus1(cases, 30)
    assert score(best, cases) > 0.03
    assert pminus1_success(e, cases) > 0.03
