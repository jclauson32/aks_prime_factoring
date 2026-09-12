"""Factoring from the AKS binomial identity, one Pascal row at a time.

Row ``n`` of Pascal's triangle reduced mod ``n`` has an interior that vanishes
identically iff ``n`` is prime.  When it does not vanish, *where* it first fails
and *what value* it fails with together name a factor of ``n``:

    the first non-zero interior entry of row ``n`` mod ``n`` sits at
    ``k = spf(n)`` and equals ``n / spf(n)``.

See ``docs/THEORY.md`` for proofs and ``docs/FINDINGS.md`` for what this does
and does not buy you computationally.
"""

from .arith import factorize, is_prime, multiplicative_order
from .central import central_split, central_sum, jacobi
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
from .cfrac import cfrac
from .divseq import elliptic_triangle_factor
from .ecm import ecm
from .batch import batch_gcd, batch_smooth
from .hyperbola import binary_search_factor, hyperbola_factor
from .implicit import implicit_factor
from .leakage import feature_accuracy, r4, tau_sum_candidates, tau_table
from .nfs import number_field_sieve
from .qs import quadratic_sieve
from .shor import shor
from .squfof import squfof
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
    "__version__",
    "aks_pow",
    "alignment_gain",
    "batch_gcd",
    "batch_smooth",
    "binary_search_factor",
    "central_split",
    "central_sum",
    "certificate",
    "cfrac",
    "check_all",
    "coppersmith_small_root",
    "ecm",
    "elliptic_triangle_factor",
    "factor",
    "factor_with_hint",
    "factorial_mod",
    "factorize",
    "fast_spf",
    "fast_split",
    "feature_accuracy",
    "first_nonzero",
    "first_zero_row",
    "fold_attack",
    "fold_coefficients",
    "fold_identity",
    "fold_norm",
    "fold_norm_expected",
    "fractal_dimension",
    "generic_rate",
    "harvey_factor",
    "harvey_search",
    "hyperbola_factor",
    "implicit_factor",
    "is_prime",
    "jacobi",
    "lehman_recover",
    "lll",
    "lucas_v",
    "multiplicative_order",
    "multipoint_eval",
    "norm_one_search",
    "number_field_sieve",
    "pascal_spf",
    "pascal_split",
    "pollard_pminus1",
    "q_first_hit",
    "q_pascal_row",
    "q_period",
    "quadratic_sieve",
    "r4",
    "render",
    "residue_shape",
    "row_entry",
    "row_exact",
    "row_prefix",
    "row_series",
    "row_support",
    "schnorr_lenstra_split",
    "shift_split",
    "shor",
    "squfof",
    "support_count",
    "tau_sum_candidates",
    "tau_table",
    "threshold_spf",
    "williams_pplus1",
    "zero_divisor_density",
]
