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
from .classgroup import schnorr_lenstra_split
from .cyclo import (
    lucas_v,
    norm_one_search,
    pollard_pminus1,
    williams_pplus1,
)
from .generic import alignment_gain, generic_rate, zero_divisor_density
from .lattice import coppersmith_small_root, factor_with_hint, lll
from .harvey import harvey_factor, harvey_search, lehman_recover
from .fractal import (
    first_zero_row,
    fractal_dimension,
    render,
    support_count,
)
from .factor import certificate, factor, pascal_split, pascal_spf
from .fast import (
    factorial_mod,
    fast_split,
    fast_spf,
    multipoint_eval,
    threshold_spf,
)
from .pascal import (
    first_nonzero,
    residue_shape,
    row_entry,
    row_exact,
    row_prefix,
    row_series,
    row_support,
)
from .qpascal import q_first_hit, q_pascal_row, q_period, shift_split
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
    "factorial_mod",
    "fast_split",
    "fast_spf",
    "first_nonzero",
    "first_zero_row",
    "fractal_dimension",
    "harvey_factor",
    "harvey_search",
    "lehman_recover",
    "lll",
    "generic_rate",
    "zero_divisor_density",
    "alignment_gain",
    "coppersmith_small_root",
    "factor_with_hint",
    "fold_attack",
    "fold_coefficients",
    "fold_identity",
    "fold_norm",
    "fold_norm_expected",
    "is_prime",
    "lucas_v",
    "norm_one_search",
    "multipoint_eval",
    "pascal_split",
    "pascal_spf",
    "pollard_pminus1",
    "q_first_hit",
    "q_pascal_row",
    "q_period",
    "williams_pplus1",
    "residue_shape",
    "row_entry",
    "row_exact",
    "row_prefix",
    "row_series",
    "render",
    "row_support",
    "support_count",
    "threshold_spf",
    "schnorr_lenstra_split",
    "shift_split",
    "__version__",
]
