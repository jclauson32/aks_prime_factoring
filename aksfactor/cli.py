"""Command line interface: ``python -m aksfactor <command> ...``."""

from __future__ import annotations

import argparse
import sys
from math import gcd

from .arith import is_prime
from .factor import certificate, factor, pascal_split
from .pascal import row_entry, row_prefix, row_support
from .classgroup import schnorr_lenstra_split
from .cyclo import norm_one_search, pollard_pminus1
from .fractal import first_zero_row, fractal_dimension, render, support_count
from .qpascal import q_first_hit, q_period, shift_split
from .ring import fold_attack
from .theorems import check_all


def _cmd_factor(args) -> int:
    n = args.n
    if n < 2:
        print("n must be >= 2")
        return 1
    if is_prime(n):
        print(f"{n} is prime -- the interior of row {n} mod {n} is identically zero.")
        return 0
    facs = factor(n, bound=args.bound, mode=args.mode)
    pretty = " * ".join(
        f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(facs.items())
    )
    print(f"{n} = {pretty}")
    split = pascal_split(n, bound=args.bound, mode=args.mode)
    if split:
        cert = certificate(n, split[0])
        print(f"  certificate: {cert['identity']}")
        print(f"  first non-zero residue of row {n} is at k = {cert['position']}")
        print(f"  smallest prime factor {cert['p']}, cofactor {cert['cofactor']}")
    return 0


def _cmd_row(args) -> int:
    n, kmax = args.n, args.kmax
    row = row_prefix(n, kmax)
    width = max(len(str(v)) for v in row)
    for k, v in enumerate(row):
        flag = ""
        if 0 < k < n and v:
            d = gcd(n, k)
            flag = f"   <- x = n/gcd(n,k) = {n // d}, y = {v // (n // d)}"
        print(f"  k={k:<6d} C(n,k) mod n = {v:>{width}d}{flag}")
    sup = [k for k in range(1, min(kmax, n - 1) + 1) if row[k]]
    print(f"\nsupport within k <= {kmax}: {sup if len(sup) <= 40 else sup[:40] + ['...']}")
    return 0


def _cmd_entry(args) -> int:
    r = row_entry(args.n, args.k)
    d = gcd(args.n, args.k)
    print(f"C({args.n}, {args.k}) mod {args.n} = {r}")
    if r:
        print(f"  x = n/gcd(n,k) = {args.n // d}, y = {r // (args.n // d)}")
        print(f"  gcd(residue, n) = {gcd(r, args.n)}  (a proper divisor of n)")
    return 0


def _cmd_verify(args) -> int:
    failures = 0
    for n in range(4, args.upto + 1):
        for name, (ok, detail) in check_all(n, kmax=args.kmax).items():
            if not ok:
                failures += 1
                print(f"FAIL n={n} {name}: {detail}")
    print(f"checked n = 4..{args.upto}: {'all theorems hold' if not failures else f'{failures} failures'}")
    return 1 if failures else 0


def _cmd_fold(args) -> int:
    res = fold_attack(args.n, rmax=args.rmax, bases=tuple(args.bases))
    if res["factor"]:
        print(f"fold attack found {res['factor']} * {res['cofactor']} = {args.n}")
        print(f"  at r={res['r']}, j={res['j']}, a={res['a']} after {res['pairs']} coefficients")
    else:
        print(f"no factor from folds r <= {args.rmax} ({res['pairs']} coefficients inspected)")
    return 0


def _cmd_normone(args) -> int:
    n = args.n
    pm1 = pollard_pminus1(n, bound=args.bound)
    found = norm_one_search(n, bound=args.bound, degree=args.degree,
                            bases=range(3, args.bases))
    print(f"n = {n}")
    print(f"  Pollard p-1  (bound {args.bound}): "
          f"{pm1 if pm1 else 'no factor'}")
    if found:
        print(f"  norm-one d={args.degree} (bound {args.bound}): "
              f"{found[0]}  [base {found[1]}]")
        print(f"  {n} = {found[0]} * {n // found[0]}")
    else:
        print(f"  norm-one d={args.degree} (bound {args.bound}): no factor")
    return 0


def _cmd_qrow(args) -> int:
    n, q = args.n, args.base
    print(f"q-Pascal row of n = {n}, base q = {q}")
    if gcd(q - 1, n) != 1:
        print(f"  (base is degenerate: gcd(q-1, n) = {gcd(q - 1, n)} splits n by itself)")
    got = shift_split(n, q, args.kmax)
    if got:
        j, f = got
        print(f"  clause-1 hit at j = {j}: ord(q) divides n - {j}")
        print(f"  factor {f}, cofactor {n // f}")
    else:
        print(f"  no clause-1 hit for j <= {args.kmax}")
    return 0


def _cmd_classgroup(args) -> int:
    n = args.n
    got = schnorr_lenstra_split(n, multipliers=range(1, args.kmax + 1),
                                bound=args.bound, forms_per_disc=args.forms)
    print(f"n = {n}")
    if got:
        print(f"  class group Cl(-{got[1]}*n) yielded an ambiguous form")
        print(f"  {n} = {got[0]} * {n // got[0]}   (multiplier k = {got[1]})")
    else:
        print(f"  no ambiguous form found for k <= {args.kmax}, bound {args.bound}")
    return 0


def _cmd_plot(args) -> int:
    n = args.n
    legend = tuple(args.legend) if args.legend else None
    for line in render(args.rows, n, legend=legend):
        print(line)
    if legend:
        print(f"\n  p={legend[0]} kills 'p',  q={legend[1]} kills 'q',  "
              f"'0' = zero mod {n} (both),  '#' = survivor")
    if is_prime(n):
        print(f"\n  mod {n}: Sierpinski gasket of ratio {n}, "
              f"box dimension {fractal_dimension(n):.4f}")
        print(f"  survivors in the first {args.rows} rows: "
              f"{support_count(args.rows, n):,}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="aksfactor", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("factor", help="factor n via the first non-zero Pascal residue")
    p.add_argument("n", type=int)
    p.add_argument("--bound", type=int, default=None)
    p.add_argument("--mode", choices=("certified", "scan", "fast"), default="certified")
    p.set_defaults(func=_cmd_factor)

    p = sub.add_parser("row", help="print a prefix of row n mod n")
    p.add_argument("n", type=int)
    p.add_argument("--kmax", type=int, default=40)
    p.set_defaults(func=_cmd_row)

    p = sub.add_parser("entry", help="one entry C(n,k) mod n, for n of any size")
    p.add_argument("n", type=int)
    p.add_argument("k", type=int)
    p.set_defaults(func=_cmd_entry)

    p = sub.add_parser("verify", help="run every theorem checker over a range")
    p.add_argument("--upto", type=int, default=400)
    p.add_argument("--kmax", type=int, default=None)
    p.set_defaults(func=_cmd_verify)

    p = sub.add_parser("fold", help="run the AKS-ring fold attack")
    p.add_argument("n", type=int)
    p.add_argument("--rmax", type=int, default=200)
    p.add_argument("--bases", type=int, nargs="*", default=[1])
    p.set_defaults(func=_cmd_fold)

    p = sub.add_parser("normone",
                       help="norm-one cyclotomic factoring vs Pollard p-1")
    p.add_argument("n", type=int)
    p.add_argument("--bound", type=int, default=5000)
    p.add_argument("--degree", type=int, default=2)
    p.add_argument("--bases", type=int, default=60)
    p.set_defaults(func=_cmd_normone)

    p = sub.add_parser("qrow", help="q-deformed Pascal row: tunable-period search")
    p.add_argument("n", type=int)
    p.add_argument("--base", type=int, default=2)
    p.add_argument("--kmax", type=int, default=10000)
    p.set_defaults(func=_cmd_qrow)

    p = sub.add_parser("classgroup",
                       help="Schnorr-Lenstra: factor via ambiguous forms in Cl(-kn)")
    p.add_argument("n", type=int)
    p.add_argument("--kmax", type=int, default=30)
    p.add_argument("--bound", type=int, default=200)
    p.add_argument("--forms", type=int, default=3)
    p.set_defaults(func=_cmd_classgroup)

    p = sub.add_parser("plot", help="render Pascal's triangle mod n as a fractal")
    p.add_argument("n", type=int)
    p.add_argument("--rows", type=int, default=27)
    p.add_argument("--legend", type=int, nargs=2, default=None,
                   metavar=("P", "Q"))
    p.set_defaults(func=_cmd_plot)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
