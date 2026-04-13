#!/usr/bin/env python3
"""Generate DD (double-double) polynomial coefficients for multifloats.

Single source of truth for all polynomial and conversion constants used
by both src/multifloats.hh (C++) and fsrc/multifloats.fypp (Fortran).

Writes:
    src/dd_constants.hh      — C++ inline constexpr arrays
    src/dd_constants.f90.inc — Fortran real(dp), parameter arrays

Usage:
    python3 scripts/gen_constants.py          # generate both files
    python3 scripts/gen_constants.py --check  # verify without writing

Requirements: mpmath (pip install mpmath)
"""

import os
import sys
from mpmath import mp, mpf, log, log10, factorial, pi
from mpmath import bernoulli, sqrt

mp.dps = 60  # 60 decimal digits — well above DD's ~32

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
CPP_OUT = os.path.join(PROJECT_DIR, "src", "dd_constants.hh")
F90_OUT = os.path.join(PROJECT_DIR, "src", "dd_constants.f90.inc")


# =============================================================================
# Utility
# =============================================================================

def to_dd(value):
    """Split an mpf value into (hi, lo) double-double pair."""
    hi = float(value)
    lo = float(value - mpf(hi))
    return hi, lo


def verify_dd(value, hi, lo, name=""):
    """Verify |value - (hi+lo)| / |value| < tol. Returns relative error."""
    dd = mpf(hi) + mpf(lo)
    if value == 0:
        return abs(float(dd))
    return float(abs((dd - value) / value))


# =============================================================================
# Constant definitions
# =============================================================================

def gen_exp2_coefs(n=14):
    """c[k] = (ln2)^k / k!"""
    ln2 = log(mpf(2))
    return [ln2**k / factorial(k) for k in range(n)]


def gen_log2_coefs(n):
    """c[k] = 2/((2k+1)*ln2)"""
    ln2 = log(mpf(2))
    return [mpf(2) / ((2*k + 1) * ln2) for k in range(n)]


def gen_log2_table(n=32):
    """log2 lookup: centers[i] = 1 + (2i+1)/64, values[i] = log2(centers[i])."""
    centers = [1.0 + (2*i + 1) / 64.0 for i in range(n)]
    values = [log(mpf(c)) / log(mpf(2)) for c in centers]
    return centers, values


def gen_sin_taylor(n=13):
    """sin(x)/x Taylor: c[k] = (-1)^k / (2k+1)!"""
    return [(-1)**k / factorial(2*k + 1) for k in range(n)]


def gen_cos_taylor(n=13):
    """cos(x) Taylor: c[k] = (-1)^k / (2k)!"""
    return [(-1)**k / factorial(2*k) for k in range(n)]


def gen_sinpi_coefs(n=15):
    """c[k] = (-1)^k * pi^(2k+1) / (2k+1)!"""
    return [(-1)**k * pi**(2*k + 1) / factorial(2*k + 1) for k in range(n)]


def gen_cospi_coefs(n=15):
    """c[k] = (-1)^k * pi^(2k) / (2k)!"""
    return [(-1)**k * pi**(2*k) / factorial(2*k) for k in range(n)]


def gen_sinh_taylor(n=9):
    """sinh(x)/x Taylor: c[k] = 1 / (2k+1)!"""
    return [mpf(1) / factorial(2*k + 1) for k in range(n)]


def gen_asinh_taylor(n=15):
    """c[k] = (-1)^k (2k)! / (4^k (k!)^2 (2k+1))"""
    return [(-1)**k * factorial(2*k) / (4**k * factorial(k)**2 * (2*k + 1))
            for k in range(n)]


def gen_atanh_taylor(n=15):
    """c[k] = 1/(2k+1)"""
    return [mpf(1) / (2*k + 1) for k in range(n)]


def gen_stirling_coefs(n=13):
    """c[k] = B_{2k} / (2k*(2k-1)), k=1..n"""
    return [bernoulli(2*k) / (2*k * (2*k - 1)) for k in range(1, n + 1)]


# =============================================================================
# Collect all constants
# =============================================================================

def collect_all():
    """Return ordered list of constant groups. Each is a dict with:
       kind: 'scalar' | 'array' | 'cw3' | 'dp_array'
       name, comment, and value data.
    """
    groups = []

    def scalar(name, value, comment):
        hi, lo = to_dd(value)
        groups.append(dict(kind='scalar', name=name, hi=hi, lo=lo,
                           exact=value, comment=comment))

    def array(name, values, comment):
        pairs = [to_dd(v) for v in values]
        groups.append(dict(kind='array', name=name,
                           hi=[p[0] for p in pairs],
                           lo=[p[1] for p in pairs],
                           exact=values, comment=comment))

    def dp_array(name, values, comment):
        """Single-precision (dp only) array — no lo part."""
        groups.append(dict(kind='dp_array', name=name,
                           values=values, comment=comment))

    # --- Conversion constants ---
    scalar('log2_e',       1 / log(mpf(2)),        'log2(e)')
    scalar('ln_2',         log(mpf(2)),            'ln(2)')
    scalar('log10_2',      log10(mpf(2)),          'log10(2)')
    scalar('inv_pi',       1 / pi,                 '1/pi')
    scalar('pi_dd',        pi,                     'pi')
    scalar('half_log_2pi', log(2 * pi) / 2,        '(1/2)*log(2*pi)')
    scalar('log_pi',       log(pi),                'log(pi)')

    # pi/2 Cody-Waite 3-part
    p = pi / 2
    cw1 = float(p)
    cw2 = float(p - mpf(cw1))
    cw3 = float(p - mpf(cw1) - mpf(cw2))
    groups.append(dict(kind='cw3', name='pi_half_cw',
                       v1=cw1, v2=cw2, v3=cw3,
                       comment='pi/2 Cody-Waite 3-part (~161 bits)'))

    # erf efx
    scalar('erf_efx', 2 / sqrt(pi) - 1, '2/sqrt(pi) - 1')
    # erf(1) — this is intentionally just a hi part (rounded to exact dp)
    groups.append(dict(kind='scalar', name='erf_const',
                       hi=0.845062911510467529296875, lo=0.0,
                       exact=mpf('0.845062911510467529296875'),
                       comment='erf(1) rounded to exact dp'))

    # --- exp2 ---
    array('exp2_coefs', gen_exp2_coefs(14), 'exp2: c[k] = (ln2)^k / k!')
    groups.append(dict(kind='exp2_clamp', comment='exp2 input clamps'))

    # --- log2 ---
    array('log2_narrow', gen_log2_coefs(7),
          'log2 narrow: c[k] = 2/((2k+1)*ln2)')
    array('log2_wide', gen_log2_coefs(9),
          'log2 wide: c[k] = 2/((2k+1)*ln2)')

    # log2 table
    centers, values = gen_log2_table(32)
    dp_array('log2_centers', centers, 'log2 lookup centers')
    array('log2_values', values, 'log2 lookup values')

    # --- trig ---
    array('sin_taylor', gen_sin_taylor(13),
          'sin(x)/x Taylor: c[k] = (-1)^k / (2k+1)!')
    array('cos_taylor', gen_cos_taylor(13),
          'cos(x) Taylor: c[k] = (-1)^k / (2k)!')
    array('sinpi_coefs', gen_sinpi_coefs(15),
          'sinpi: c[k] = (-1)^k * pi^(2k+1) / (2k+1)!')
    array('cospi_coefs', gen_cospi_coefs(15),
          'cospi: c[k] = (-1)^k * pi^(2k) / (2k)!')

    # --- hyperbolic ---
    array('sinh_taylor', gen_sinh_taylor(9),
          'sinh(x)/x Taylor: c[k] = 1/(2k+1)!')
    array('asinh_taylor', gen_asinh_taylor(15),
          'asinh(x)/x Taylor: c[k] = (-1)^k(2k)!/(4^k(k!)^2(2k+1))')
    array('atanh_taylor', gen_atanh_taylor(15),
          'atanh(x)/x Taylor: c[k] = 1/(2k+1)')

    # --- gamma ---
    array('stirling_coefs', gen_stirling_coefs(13),
          'Stirling: c[k] = B_{2k}/(2k*(2k-1)), k=1..13')

    return groups


# =============================================================================
# Output: C++
# =============================================================================

def write_cpp(groups, f):
    f.write("// ==========================================================================\n")
    f.write("// AUTO-GENERATED by scripts/gen_constants.py — do not edit manually.\n")
    f.write("// Regenerate: python3 scripts/gen_constants.py\n")
    f.write("// ==========================================================================\n")
    f.write("#pragma once\n\n")

    for g in groups:
        kind = g['kind']
        if kind == 'scalar':
            f.write(f"// {g['comment']}\n")
            f.write(f"inline constexpr double {g['name']}_hi = {g['hi']:23.17e};\n")
            f.write(f"inline constexpr double {g['name']}_lo = {g['lo']:23.17e};\n\n")

        elif kind == 'cw3':
            f.write(f"// {g['comment']}\n")
            f.write(f"inline constexpr double {g['name']}1 = {g['v1']:23.17e};\n")
            f.write(f"inline constexpr double {g['name']}2 = {g['v2']:23.17e};\n")
            f.write(f"inline constexpr double {g['name']}3 = {g['v3']:23.17e};\n\n")

        elif kind == 'exp2_clamp':
            f.write("// exp2 input clamps\n")
            f.write("inline constexpr double exp2_min_d = -1022.0;\n")
            f.write("inline constexpr double exp2_max_d = 1023.9999999999998;\n\n")

        elif kind == 'array':
            n = len(g['hi'])
            f.write(f"// {g['comment']}\n")
            f.write(f"inline constexpr double {g['name']}_hi[{n}] = {{\n")
            for i in range(0, n, 2):
                vals = [f"{g['hi'][j]:23.17e}" for j in range(i, min(i+2, n))]
                end = "};" if i + 2 >= n else ","
                f.write(f"    {', '.join(vals)}{end}\n")
            f.write(f"inline constexpr double {g['name']}_lo[{n}] = {{\n")
            for i in range(0, n, 2):
                vals = [f"{g['lo'][j]:23.17e}" for j in range(i, min(i+2, n))]
                end = "};" if i + 2 >= n else ","
                f.write(f"    {', '.join(vals)}{end}\n")
            f.write("\n")

        elif kind == 'dp_array':
            vals = g['values']
            n = len(vals)
            f.write(f"// {g['comment']}\n")
            f.write(f"inline constexpr double {g['name']}[{n}] = {{\n")
            for i in range(0, n, 4):
                chunk = [str(v) for v in vals[i:i+4]]
                end = "};" if i + 4 >= n else ","
                f.write(f"    {', '.join(chunk)}{end}\n")
            f.write("\n")


# =============================================================================
# Output: Fortran
# =============================================================================

def write_f90(groups, f):
    f.write("    ! ======================================================================\n")
    f.write("    ! AUTO-GENERATED by scripts/gen_constants.py — do not edit manually.\n")
    f.write("    ! Regenerate: python3 scripts/gen_constants.py\n")
    f.write("    ! ======================================================================\n\n")

    for g in groups:
        kind = g['kind']
        if kind == 'scalar':
            f.write(f"    ! {g['comment']}\n")
            f.write(f"    real(dp), parameter :: {g['name']}_hi = {g['hi']:23.17e}_dp\n")
            f.write(f"    real(dp), parameter :: {g['name']}_lo = {g['lo']:23.17e}_dp\n\n")

        elif kind == 'cw3':
            f.write(f"    ! {g['comment']}\n")
            f.write(f"    real(dp), parameter :: {g['name']}1 = {g['v1']:23.17e}_dp\n")
            f.write(f"    real(dp), parameter :: {g['name']}2 = {g['v2']:23.17e}_dp\n")
            f.write(f"    real(dp), parameter :: {g['name']}3 = {g['v3']:23.17e}_dp\n\n")

        elif kind == 'exp2_clamp':
            f.write("    ! exp2 input clamps\n")
            f.write("    real(dp), parameter :: exp2_min = -1022.0_dp\n")
            f.write("    real(dp), parameter :: exp2_max = 1023.9999999999998_dp\n\n")

        elif kind == 'array':
            n = len(g['hi'])
            f.write(f"    ! {g['comment']}\n")
            f.write(f"    real(dp), parameter :: {g['name']}_hi({n}) = [ &\n")
            for i in range(0, n, 2):
                vals = [f"{g['hi'][j]:23.17e}_dp" for j in range(i, min(i+2, n))]
                end = " ]" if i + 2 >= n else ", &"
                f.write(f"        {', '.join(vals)}{end}\n")
            f.write(f"    real(dp), parameter :: {g['name']}_lo({n}) = [ &\n")
            for i in range(0, n, 2):
                vals = [f"{g['lo'][j]:23.17e}_dp" for j in range(i, min(i+2, n))]
                end = " ]" if i + 2 >= n else ", &"
                f.write(f"        {', '.join(vals)}{end}\n")
            f.write("\n")

        elif kind == 'dp_array':
            vals = g['values']
            n = len(vals)
            f.write(f"    ! {g['comment']}\n")
            f.write(f"    real(dp), parameter :: {g['name']}({n}) = [ &\n")
            for i in range(0, n, 4):
                chunk = [f"{v}_dp" for v in vals[i:i+4]]
                end = " ]" if i + 4 >= n else ", &"
                f.write(f"        {', '.join(chunk)}{end}\n")
            f.write("\n")


# =============================================================================
# Verification
# =============================================================================

def verify_all(groups):
    max_err = 0.0
    n_checked = 0
    for g in groups:
        if g['kind'] == 'scalar':
            err = verify_dd(g['exact'], g['hi'], g['lo'], g['name'])
            max_err = max(max_err, err)
            n_checked += 1
        elif g['kind'] == 'array':
            for i in range(len(g['hi'])):
                err = verify_dd(g['exact'][i], g['hi'][i], g['lo'][i],
                                f"{g['name']}[{i}]")
                max_err = max(max_err, err)
                n_checked += 1
    return max_err, n_checked


# =============================================================================
# Main
# =============================================================================

def main():
    check_only = '--check' in sys.argv

    groups = collect_all()
    max_err, n_checked = verify_all(groups)

    print(f"Verified {n_checked} constants, max DD conversion error: {max_err:.3e}",
          file=sys.stderr)

    if max_err > 1e-31:
        print("WARNING: some constants exceed 1e-31 relative error", file=sys.stderr)

    if check_only:
        print("Check-only mode, not writing files.", file=sys.stderr)
        return

    with open(CPP_OUT, 'w') as f:
        write_cpp(groups, f)
    print(f"Wrote {CPP_OUT}", file=sys.stderr)

    with open(F90_OUT, 'w') as f:
        write_f90(groups, f)
    print(f"Wrote {F90_OUT}", file=sys.stderr)


if __name__ == '__main__':
    main()
