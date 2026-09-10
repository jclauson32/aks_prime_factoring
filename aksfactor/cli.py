"""Command line interface: ``python -m aksfactor <command> ...``."""

from __future__ import annotations

import argparse
import sys
from math import gcd

from .arith import is_prime
from .factor import certificate, factor, pascal_split
from .pascal import row_entry, row_prefix, row_support
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

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
