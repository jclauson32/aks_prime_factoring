"""Exhaustive verification of T1-T5 and the original conjecture."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aksfactor.theorems import (
    check_all,
    check_gist,
    check_t1_shape,
    check_t2_coprime_vanishing,
    check_t3_exact_value,
    check_t4_first_nonzero,
    check_t5_closed_form,
)

RANGE = range(4, 320)


def test_t1_shape():
    for n in RANGE:
        ok, detail = check_t1_shape(n)
        assert ok, detail


def test_t2_coprime_vanishing():
    for n in RANGE:
        ok, detail = check_t2_coprime_vanishing(n)
        assert ok, detail


def test_t3_exact_value():
    for n in RANGE:
        ok, detail = check_t3_exact_value(n)
        assert ok, detail


def test_t4_first_nonzero():
    for n in RANGE:
        ok, detail = check_t4_first_nonzero(n)
        assert ok, detail


def test_t5_closed_form():
    for n in range(4, 160):
        ok, detail = check_t5_closed_form(n)
        assert ok, detail


def test_gist_holds():
    for n in RANGE:
        ok, detail = check_gist(n)
        assert ok, detail


def test_check_all_bundles():
    for n in (12, 35, 97, 561, 1024, 2047):
        for name, (ok, detail) in check_all(n, kmax=200).items():
            assert ok, (n, name, detail)
