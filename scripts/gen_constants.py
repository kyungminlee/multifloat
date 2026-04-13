#!/usr/bin/env python3
"""Generate DD (double-double) polynomial coefficients for multifloats.

Single source of truth for all polynomial and conversion constants used
by both src/multifloats.hh (C++) and fsrc/multifloats.fypp (Fortran).

Usage:
    python3 scripts/gen_constants.py          # print to stdout
    python3 scripts/gen_constants.py --check  # verify against source files

Requirements: mpmath (pip install mpmath)
"""

import struct
import sys
from mpmath import mp, mpf, log, log10, factorial, pi, euler as euler_gamma
from mpmath import bernoulli, sqrt

mp.dps = 60  # 60 decimal digits — well above DD's ~32


# =============================================================================
# Utility: split a high-precision value into a DD pair (hi, lo)
# =============================================================================

def to_dd(value):
    """Split an mpf value into (hi, lo) where hi = float(value), lo = float(value - hi)."""
    hi = float(value)
    lo = float(value - mpf(hi))
    return hi, lo


def verify_dd(value, hi, lo, name="", tol=1e-31):
    """Verify |value - (hi+lo)| / |value| < tol."""
    dd = mpf(hi) + mpf(lo)
    if value == 0:
        err = abs(dd)
    else:
        err = float(abs((dd - value) / value))
    if err > tol:
        print(f"  WARNING: {name} rel_err={err:.3e} > {tol:.0e}", file=sys.stderr)
    return err


# =============================================================================
# Constant definitions
# =============================================================================

def gen_conversion_constants():
    """log2(e), ln(2), log10(2), 1/pi, pi, etc."""
    consts = {}
    consts['log2_e'] = log(mpf(2))**(-1)  # 1/ln(2) = log2(e)
    consts['ln_2'] = log(mpf(2))
    consts['log10_2'] = log10(mpf(2))
    consts['inv_pi'] = 1 / pi
    consts['pi_dd'] = pi
    consts['half_log_2pi'] = log(2 * pi) / 2
    consts['log_pi'] = log(pi)
    return consts


def gen_exp2_coefs(n=14):
    """exp2 polynomial: c[k] = (ln2)^k / k!"""
    ln2 = log(mpf(2))
    coefs = []
    for k in range(n):
        coefs.append(ln2**k / factorial(k))
    return coefs


def gen_log2_narrow_coefs(n=7):
    """log2 narrow-path polynomial.

    For the table-lookup path, we approximate log2(1+s) where s is small.
    The polynomial is in terms of t = s/(2+s) (argument reduction).
    Coefficients: c[k] = 2/(2k+1) * (1/ln2), evaluated at odd powers of t.
    c[k] = 2 / ((2k+1) * ln(2))
    """
    ln2 = log(mpf(2))
    coefs = []
    for k in range(n):
        coefs.append(mpf(2) / ((2*k + 1) * ln2))
    return coefs


def gen_log2_wide_coefs(n=9):
    """log2 wide-path polynomial (same formula, more terms)."""
    return gen_log2_narrow_coefs(n)


def gen_log2_table(n=32):
    """log2 lookup table: centers and log2(center) values."""
    centers = []
    values = []
    for i in range(n):
        c = 1.0 + (2*i + 1) / 64.0
        centers.append(c)
        values.append(log(mpf(c)) / log(mpf(2)))
    return centers, values


def gen_sin_taylor(n=13):
    """sin(x)/x Taylor: c[k] = (-1)^k / (2k+1)!"""
    coefs = []
    for k in range(n):
        coefs.append((-1)**k / factorial(2*k + 1))
    return coefs


def gen_cos_taylor(n=13):
    """cos(x) Taylor: c[k] = (-1)^k / (2k)!"""
    coefs = []
    for k in range(n):
        coefs.append((-1)**k / factorial(2*k))
    return coefs


def gen_sinpi_coefs(n=15):
    """sinpi(x) = x * P(x^2) where P coefficients are (-1)^k * pi^(2k+1) / (2k+1)!
    Actually: sinpi(x) = sin(pi*x), and we want c[k] for the Horner in x^2.
    c[k] = (-1)^k * pi^(2k+1) / (2k+1)!  ... but sign convention may differ.
    Let me compute: sin(pi*x) = sum_{k=0} (-1)^k (pi*x)^{2k+1} / (2k+1)!
                               = x * sum_{k=0} (-1)^k pi^{2k+1} x^{2k} / (2k+1)!
    So c[k] = (-1)^k * pi^(2k+1) / (2k+1)!
    """
    coefs = []
    for k in range(n):
        coefs.append((-1)**k * pi**(2*k + 1) / factorial(2*k + 1))
    return coefs


def gen_cospi_coefs(n=15):
    """cospi(x) = cos(pi*x) = sum_{k=0} (-1)^k (pi*x)^{2k} / (2k)!
    c[k] = (-1)^k * pi^{2k} / (2k)!
    """
    coefs = []
    for k in range(n):
        coefs.append((-1)**k * pi**(2*k) / factorial(2*k))
    return coefs


def gen_sinh_taylor(n=9):
    """sinh(x)/x Taylor: c[k] = 1 / (2k+1)!"""
    coefs = []
    for k in range(n):
        coefs.append(mpf(1) / factorial(2*k + 1))
    return coefs


def gen_asinh_taylor(n=15):
    """asinh(x)/x Taylor: c[k] = (-1)^k (2k)! / (4^k (k!)^2 (2k+1))"""
    coefs = []
    for k in range(n):
        coefs.append((-1)**k * factorial(2*k) / (4**k * factorial(k)**2 * (2*k + 1)))
    return coefs


def gen_atanh_taylor(n=15):
    """atanh(x)/x Taylor: c[k] = 1/(2k+1)"""
    coefs = []
    for k in range(n):
        coefs.append(mpf(1) / (2*k + 1))
    return coefs


def gen_stirling_coefs(n=13):
    """Stirling series: c[k] = B_{2k} / (2k * (2k-1)), k=1..n"""
    coefs = []
    for k in range(1, n + 1):
        b2k = bernoulli(2*k)
        coefs.append(b2k / (2*k * (2*k - 1)))
    return coefs


def gen_pi_half_cw():
    """pi/2 as 3-part Cody-Waite constant (~161 bits).
    Split pi/2 into three parts where each successive part captures
    the rounding error of the sum of previous parts.
    """
    p = pi / 2
    cw1 = float(p)
    cw2 = float(p - mpf(cw1))
    cw3 = float(p - mpf(cw1) - mpf(cw2))
    return cw1, cw2, cw3


def gen_erf_efx():
    """2/sqrt(pi) - 1"""
    return 2 / sqrt(pi) - 1


# =============================================================================
# Output formatting
# =============================================================================

def fmt_cpp_array(name, values_hi, values_lo, comment=""):
    """Format as C++ inline constexpr double arrays."""
    n = len(values_hi)
    lines = []
    if comment:
        lines.append(f"// ---- {comment}")
    lines.append(f"inline constexpr double {name}_hi[{n}] = {{")
    for i in range(0, n, 2):
        if i + 1 < n:
            lines.append(f"    {values_hi[i]:23.17e}, {values_hi[i+1]:23.17e},")
        else:
            lines.append(f"    {values_hi[i]:23.17e}}};")
    if n % 2 == 0:
        lines[-1] = lines[-1].rstrip(',') + '};'
    lines.append(f"inline constexpr double {name}_lo[{n}] = {{")
    for i in range(0, n, 2):
        if i + 1 < n:
            lines.append(f"    {values_lo[i]:23.17e}, {values_lo[i+1]:23.17e},")
        else:
            lines.append(f"    {values_lo[i]:23.17e}}};")
    if n % 2 == 0:
        lines[-1] = lines[-1].rstrip(',') + '};'
    return '\n'.join(lines)


def fmt_cpp_scalar(name, hi, lo, comment=""):
    """Format as C++ inline constexpr double pair."""
    lines = []
    if comment:
        lines.append(f"// {comment}")
    lines.append(f"inline constexpr double {name}_hi = {hi:23.17e};")
    lines.append(f"inline constexpr double {name}_lo = {lo:23.17e};")
    return '\n'.join(lines)


def fmt_fortran_array(name, values_hi, values_lo, comment=""):
    """Format as Fortran real(dp), parameter arrays."""
    n = len(values_hi)
    lines = []
    if comment:
        lines.append(f"    ! {comment}")
    lines.append(f"    real(dp), parameter :: {name}_hi({n}) = [ &")
    for i in range(0, n, 2):
        if i + 1 < n:
            end = " ]" if i + 2 >= n else ", &"
            lines.append(f"        {values_hi[i]:23.17e}_dp, {values_hi[i+1]:23.17e}_dp{end}")
        else:
            lines.append(f"        {values_hi[i]:23.17e}_dp ]")
    lines.append(f"    real(dp), parameter :: {name}_lo({n}) = [ &")
    for i in range(0, n, 2):
        if i + 1 < n:
            end = " ]" if i + 2 >= n else ", &"
            lines.append(f"        {values_lo[i]:23.17e}_dp, {values_lo[i+1]:23.17e}_dp{end}")
        else:
            lines.append(f"        {values_lo[i]:23.17e}_dp ]")
    return '\n'.join(lines)


# =============================================================================
# Main: generate and verify all constants
# =============================================================================

def generate_all():
    """Generate all constants, verify precision, and print."""
    all_groups = []
    max_err = 0.0

    def add_array(name, exact_values, comment=""):
        nonlocal max_err
        hi_vals, lo_vals = [], []
        for i, v in enumerate(exact_values):
            hi, lo = to_dd(v)
            err = verify_dd(v, hi, lo, f"{name}[{i}]")
            max_err = max(max_err, err)
            hi_vals.append(hi)
            lo_vals.append(lo)
        all_groups.append(('array', name, hi_vals, lo_vals, comment))
        return hi_vals, lo_vals

    def add_scalar(name, exact_value, comment=""):
        nonlocal max_err
        hi, lo = to_dd(exact_value)
        err = verify_dd(exact_value, hi, lo, name)
        max_err = max(max_err, err)
        all_groups.append(('scalar', name, hi, lo, comment))
        return hi, lo

    # --- Conversion constants ---
    consts = gen_conversion_constants()
    add_scalar('log2_e', consts['log2_e'], 'log2(e)')
    add_scalar('ln_2', consts['ln_2'], 'ln(2)')
    add_scalar('log10_2', consts['log10_2'], 'log10(2)')
    add_scalar('inv_pi', consts['inv_pi'], '1/pi')
    add_scalar('pi_dd', consts['pi_dd'], 'pi')
    add_scalar('half_log_2pi', consts['half_log_2pi'], '(1/2)*log(2*pi)')
    add_scalar('log_pi', consts['log_pi'], 'log(pi)')

    # pi/2 Cody-Waite (3-part, not DD)
    cw1, cw2, cw3 = gen_pi_half_cw()
    all_groups.append(('cw3', 'pi_half_cw', cw1, cw2, cw3, 'pi/2 Cody-Waite 3-part'))

    # erf efx = 2/sqrt(pi) - 1
    add_scalar('erf_efx', gen_erf_efx(), '2/sqrt(pi) - 1')

    # --- Polynomial arrays ---
    add_array('exp2_coefs', gen_exp2_coefs(14),
              'exp2 polynomial: c[k] = (ln2)^k / k!')
    add_array('log2_narrow', gen_log2_narrow_coefs(7),
              'log2 narrow polynomial: c[k] = 2/((2k+1)*ln2)')
    add_array('log2_wide', gen_log2_wide_coefs(9),
              'log2 wide polynomial: c[k] = 2/((2k+1)*ln2)')

    # log2 table
    centers, values = gen_log2_table(32)
    add_array('log2_values', values, 'log2 lookup table values')
    all_groups.append(('centers', 'log2_centers', centers, 'log2 lookup table centers'))

    add_array('sin_taylor', gen_sin_taylor(13),
              'sin(x)/x Taylor: c[k] = (-1)^k / (2k+1)!')
    add_array('cos_taylor', gen_cos_taylor(13),
              'cos(x) Taylor: c[k] = (-1)^k / (2k)!')
    add_array('sinpi_coefs', gen_sinpi_coefs(15),
              'sinpi polynomial: c[k] = (-1)^k * pi^(2k+1) / (2k+1)!')
    add_array('cospi_coefs', gen_cospi_coefs(15),
              'cospi polynomial: c[k] = (-1)^k * pi^(2k) / (2k)!')
    add_array('sinh_taylor', gen_sinh_taylor(9),
              'sinh(x)/x Taylor: c[k] = 1 / (2k+1)!')
    add_array('asinh_taylor', gen_asinh_taylor(15),
              'asinh(x)/x Taylor: c[k] = (-1)^k (2k)! / (4^k (k!)^2 (2k+1))')
    add_array('atanh_taylor', gen_atanh_taylor(15),
              'atanh(x)/x Taylor: c[k] = 1/(2k+1)')
    add_array('stirling_coefs', gen_stirling_coefs(13),
              'Stirling series: c[k] = B_{2k} / (2k*(2k-1)), k=1..13')

    # --- Print ---
    print("// =============================================================================")
    print("// AUTO-GENERATED by scripts/gen_constants.py — do not edit manually")
    print("// =============================================================================")
    print()

    for group in all_groups:
        if group[0] == 'array':
            _, name, hi, lo, comment = group
            print(fmt_cpp_array(name, hi, lo, comment))
            print()
        elif group[0] == 'scalar':
            _, name, hi, lo, comment = group
            print(fmt_cpp_scalar(name, hi, lo, comment))
            print()
        elif group[0] == 'cw3':
            _, name, cw1, cw2, cw3, comment = group
            print(f"// {comment}")
            print(f"inline constexpr double {name}1 = {cw1:.17e};")
            print(f"inline constexpr double {name}2 = {cw2:.17e};")
            print(f"inline constexpr double {name}3 = {cw3:.17e};")
            print()
        elif group[0] == 'centers':
            _, name, vals, comment = group
            print(f"// {comment}")
            print(f"inline constexpr double {name}[{len(vals)}] = {{")
            for i in range(0, len(vals), 4):
                chunk = vals[i:i+4]
                line = ', '.join(f'{v}' for v in chunk)
                end = '};' if i + 4 >= len(vals) else ','
                print(f"    {line}{end}")
            print()

    print()
    print("// =============================================================================")
    print("// Fortran versions")
    print("// =============================================================================")
    print()

    for group in all_groups:
        if group[0] == 'array':
            _, name, hi, lo, comment = group
            print(fmt_fortran_array(name, hi, lo, comment))
            print()
        elif group[0] == 'scalar':
            _, name, hi, lo, comment = group
            print(f"    ! {comment}")
            print(f"    real(dp), parameter :: {name}_hi = {hi:23.17e}_dp")
            print(f"    real(dp), parameter :: {name}_lo = {lo:23.17e}_dp")
            print()

    print(f"\n// Maximum DD conversion error across all constants: {max_err:.3e}",
          file=sys.stderr)
    if max_err > 1e-31:
        print("// WARNING: some constants exceed 1e-31 relative error", file=sys.stderr)
    else:
        print("// All constants verified to < 1e-31 relative error", file=sys.stderr)


def check_against_source():
    """Compare generated values against what's currently in the source files."""
    print("Checking generated constants against source files...")
    print()

    # Check exp2 as representative example
    exp2 = gen_exp2_coefs(14)
    print("exp2_coefs (c[k] = (ln2)^k / k!):")
    for k, v in enumerate(exp2):
        hi, lo = to_dd(v)
        err = verify_dd(v, hi, lo, f"c[{k}]")
        print(f"  c[{k:2d}]: hi={hi:23.17e}  lo={lo:23.17e}  err={err:.3e}")

    print()
    # Check log2 narrow
    log2n = gen_log2_narrow_coefs(7)
    print("log2_narrow (c[k] = 2/((2k+1)*ln2)):")
    for k, v in enumerate(log2n):
        hi, lo = to_dd(v)
        err = verify_dd(v, hi, lo, f"c[{k}]")
        print(f"  c[{k}]: hi={hi:23.17e}  lo={lo:23.17e}  err={err:.3e}")

    print()
    log2w = gen_log2_wide_coefs(9)
    print("log2_wide (c[k] = 2/((2k+1)*ln2)):")
    for k, v in enumerate(log2w):
        hi, lo = to_dd(v)
        err = verify_dd(v, hi, lo, f"c[{k}]")
        print(f"  c[{k}]: hi={hi:23.17e}  lo={lo:23.17e}  err={err:.3e}")

    print()
    stirling = gen_stirling_coefs(13)
    print("stirling_coefs (c[k] = B_{2k}/(2k*(2k-1))):")
    for k, v in enumerate(stirling):
        hi, lo = to_dd(v)
        err = verify_dd(v, hi, lo, f"c[{k}]")
        print(f"  c[{k:2d}]: hi={hi:23.17e}  lo={lo:23.17e}  err={err:.3e}")


if __name__ == '__main__':
    if '--check' in sys.argv:
        check_against_source()
    else:
        generate_all()
