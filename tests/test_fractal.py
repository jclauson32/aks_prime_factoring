"""Pascal mod n as a fractal: Lucas geometry, box counts, and the row/column duals."""

import sys
from math import log
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.arith import sieve, spf_trial
from aksfactor.fractal import (
    digits,
    first_zero_row,
    first_zero_row_predicted,
    fractal_dimension,
    render,
    row_has_zero,
    row_nonzero_count,
    support_count,
)


def _brute_rows(rows, m):
    out, row = [], [1]
    for i in range(rows):
        if i:
            row = [1] + [(row[j - 1] + row[j]) % m for j in range(1, i)] + [1]
        out.append(row[:])
    return out


def test_digits():
    assert digits(0, 3) == [0]
    assert digits(5, 3) == [2, 1]
    assert digits(27, 3) == [0, 0, 0, 1]


def test_row_nonzero_count_matches_lucas():
    for p in (2, 3, 5, 7):
        for i, row in enumerate(_brute_rows(120, p)):
            assert row_nonzero_count(i, p) == sum(1 for v in row if v), (p, i)


def test_row_has_zero():
    for p in (2, 3, 5, 7):
        for i, row in enumerate(_brute_rows(120, p)):
            assert row_has_zero(i, p) == any(v == 0 for v in row), (p, i)


def test_support_count_digit_dp():
    for p in (2, 3, 5, 7):
        total = 0
        for i, row in enumerate(_brute_rows(200, p)):
            assert support_count(i, p) == total, (p, i)
            total += sum(1 for v in row if v)


def test_self_similarity():
    """The defining relation of the gasket."""
    for p in (2, 3, 5, 7, 11):
        for k in (1, 2, 3):
            assert support_count(p**k, p) == (p * (p + 1) // 2) ** k


def test_fractal_dimension_increases_with_p():
    dims = [fractal_dimension(p) for p in (2, 3, 5, 7, 11, 101)]
    assert dims == sorted(dims)
    assert all(1 < d < 2 for d in dims)
    assert abs(fractal_dimension(2) - log(3) / log(2)) < 1e-12


def test_first_zero_row_matches_prediction():
    primes = sieve(50)
    checked = 0
    for a in range(len(primes)):
        for b in range(a + 1, len(primes)):
            p, q = primes[a], primes[b]
            n = p * q
            if n > 1200:
                continue
            checked += 1
            assert first_zero_row(n, 6 * q + 60) == first_zero_row_predicted(p, q)
    assert checked > 50


def test_first_zero_row_is_larger_prime_iff_digit_condition():
    primes = sieve(50)
    for a in range(len(primes)):
        for b in range(a + 1, len(primes)):
            p, q = primes[a], primes[b]
            n = p * q
            if n > 1200:
                continue
            got = first_zero_row(n, 6 * q + 60)
            assert (got == q) == row_has_zero(q, p), (n, got, q)


def test_no_zero_below_larger_prime():
    """No row below q can contain a zero mod n -- one base-q digit means no holes."""
    primes = sieve(40)
    for a in range(len(primes)):
        for b in range(a + 1, len(primes)):
            p, q = primes[a], primes[b]
            n = p * q
            if n > 900:
                continue
            got = first_zero_row(n, 6 * q + 60)
            assert got is None or got >= q, (n, got, q)


def test_render_shapes():
    pic = render(9, 3)
    assert len(pic) == 9
    assert pic[-1].strip() == "#" * 9
    legend = render(9, 15, legend=(3, 5))
    assert all(ch in " #pq0" for line in legend for ch in line)


def test_column_divides_matches_direct_binomial():
    """Lucas, read down a column instead of along a row."""
    from math import comb

    from aksfactor.fractal import column_divides

    for p in (2, 3, 5, 7, 11):
        for c in range(0, 12):
            for i in range(c, 90):
                assert column_divides(i, c, p) == (comb(i, c) % p == 0), (p, c, i)


def test_column_zero_pattern_is_periodic():
    """For c < p the pattern is periodic with period p -- the sideways gasket."""
    from aksfactor.fractal import column_divides

    for p in (7, 11, 13):
        for c in range(1, p):
            pattern = [column_divides(i, c, p) for i in range(c, c + 5 * p)]
            for k in range(len(pattern) - p):
                assert pattern[k] == pattern[k + p], (p, c, k)
