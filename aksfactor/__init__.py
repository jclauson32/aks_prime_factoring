"""Factoring from the AKS binomial identity, one Pascal row at a time.

Row ``n`` of Pascal's triangle reduced mod ``n`` has an interior that vanishes
identically iff ``n`` is prime.  When it does not vanish, *where* it first fails
and *what value* it fails with together name a factor of ``n``:

    the first non-zero interior entry of row ``n`` mod ``n`` sits at
    ``k = spf(n)`` and equals ``n / spf(n)``.

See ``docs/THEORY.md`` for proofs and ``docs/FINDINGS.md`` for what this does
and does not buy you computationally.
"""

from .arith import factorize, is_prime
from .factor import certificate, factor, pascal_split, pascal_spf
from .fast import fast_split, fast_spf, multipoint_eval
from .pascal import (
    first_nonzero,
    residue_shape,
    row_entry,
    row_exact,
    row_prefix,
    row_series,
    row_support,
)
from .ring import (
    aks_pow,
    fold_attack,
    fold_coefficients,
    fold_identity,
    fold_norm,
    fold_norm_expected,
)
from .theorems import check_all

__version__ = "0.1.0"

__all__ = [
    "aks_pow",
    "certificate",
    "check_all",
    "factor",
    "factorize",
    "fast_split",
    "fast_spf",
    "first_nonzero",
    "fold_attack",
    "fold_coefficients",
    "fold_identity",
    "fold_norm",
    "fold_norm_expected",
    "is_prime",
    "multipoint_eval",
    "pascal_split",
    "pascal_spf",
    "residue_shape",
    "row_entry",
    "row_exact",
    "row_prefix",
    "row_series",
    "row_support",
    "__version__",
]
