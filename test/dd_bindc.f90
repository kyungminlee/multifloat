! DD arithmetic kernels with bind(c) ABI — same algorithm as the
! fypp-generated multifloats kernels, but using the C calling convention
! (pass/return by value in registers on ARM64).
module dd_bindc
  use, intrinsic :: iso_c_binding,   only: c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_fma
  implicit none
  private

  integer, parameter :: dp = 8

  type, bind(c), public :: dd_c
    real(c_double) :: hi, lo
  end type dd_c

  public :: fc_add, fc_sub, fc_mul, fc_div, fc_sqrt

contains

  pure function fc_add(a, b) result(res) bind(c)
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    real(dp) :: s, e, t, aprime, bprime, aerr, berr, c, d
    s = a%hi + b%hi
    if (.not. ieee_is_finite(s)) then
      res%hi = s; res%lo = 0.0_dp; return
    end if
    aprime = s - b%hi; bprime = s - aprime
    aerr = a%hi - aprime; berr = b%hi - bprime
    e = aerr + berr
    c = a%lo + b%lo
    aprime = c - b%lo; bprime = c - aprime
    aerr = a%lo - aprime; berr = b%lo - bprime
    d = aerr + berr
    t = s + c; bprime = t - s; c = c - bprime
    e = e + d + c
    res%hi = t + e; bprime = res%hi - t; res%lo = e - bprime
  end function

  pure function fc_sub(a, b) result(res) bind(c)
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    type(dd_c) :: nb
    nb%hi = -b%hi; nb%lo = -b%lo
    res = fc_add(a, nb)
  end function

  pure function fc_mul(a, b) result(res) bind(c)
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    real(dp) :: p00, e00, p01, p10, bprime
    p00 = a%hi * b%hi
    if (.not. ieee_is_finite(p00)) then
      res%hi = p00; res%lo = 0.0_dp; return
    end if
    e00 = ieee_fma(a%hi, b%hi, -p00)
    p01 = a%hi * b%lo; p10 = a%lo * b%hi
    e00 = e00 + (p01 + p10)
    res%hi = p00 + e00; bprime = res%hi - p00; res%lo = e00 - bprime
  end function

  pure function fc_div(a, b) result(res) bind(c)
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    real(dp) :: s, ux, p00, e00, bprime, q0, q1l
    real(dp) :: qy0, qy1, r0, r1, aprime, aerr, berr, c, t, e
    s = a%hi / b%hi
    if (.not. ieee_is_finite(s)) then
      res%hi = s; res%lo = 0.0_dp; return
    end if
    ux = 1.0_dp / b%hi
    p00 = a%hi * ux
    e00 = ieee_fma(a%hi, ux, -p00) + a%lo * ux
    q0 = p00 + e00; bprime = q0 - p00; q1l = e00 - bprime
    p00 = q0 * b%hi
    e00 = ieee_fma(q0, b%hi, -p00) + q0 * b%lo + q1l * b%hi
    qy0 = p00 + e00; bprime = qy0 - p00; qy1 = e00 - bprime
    s = qy0 + (-a%hi)
    aprime = s - (-a%hi); bprime = s - aprime
    aerr = qy0 - aprime; berr = (-a%hi) - bprime; e = aerr + berr
    c = qy1 + (-a%lo); t = s + c; bprime = t - s; c = c - bprime; e = e + c
    r0 = t + e; bprime = r0 - t; r1 = e - bprime
    p00 = r0 * ux; e00 = ieee_fma(r0, ux, -p00) + r1 * ux
    s = q0 + (-p00)
    aprime = s - (-p00); bprime = s - aprime
    aerr = q0 - aprime; berr = (-p00) - bprime; e = aerr + berr
    c = q1l + (-e00); t = s + c; bprime = t - s; c = c - bprime; e = e + c
    res%hi = t + e; bprime = res%hi - t; res%lo = e - bprime
  end function

  pure function fc_sqrt(a) result(res) bind(c)
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    real(dp) :: x0, x1, leading, u, half_u
    real(dp) :: r0, r1, p00, e00, bprime, s, c, t, e, aprime, aerr, berr
    x0 = a%hi; x1 = a%lo
    leading = sqrt(x0)
    if (.not. (x0 > 0.0_dp .and. ieee_is_finite(leading))) then
      res%hi = leading; res%lo = 0.0_dp; return
    end if
    u = 1.0_dp / leading; half_u = 0.5_dp * u
    p00 = x0 * u; e00 = ieee_fma(x0, u, -p00) + x1 * u
    r0 = p00 + e00; bprime = r0 - p00; r1 = e00 - bprime
    p00 = r0 * r0; e00 = ieee_fma(r0, r0, -p00) + 2.0_dp * (r0 * r1)
    s = p00 + e00; bprime = s - p00; c = e00 - bprime
    t = s + (-x0); aprime = t - (-x0); bprime = t - aprime
    aerr = s - aprime; berr = (-x0) - bprime; e = aerr + berr; e = e + (c - x1)
    p00 = t + e; bprime = p00 - t; e00 = e - bprime
    s = p00 * half_u; c = ieee_fma(p00, half_u, -s) + e00 * half_u
    p00 = s + c; bprime = p00 - s; e00 = c - bprime
    s = r0 + (-p00); aprime = s - (-p00); bprime = s - aprime
    aerr = r0 - aprime; berr = (-p00) - bprime; e = aerr + berr
    c = r1 + (-e00); t = s + c; bprime = t - s; c = c - bprime; e = e + c
    res%hi = t + e; bprime = res%hi - t; res%lo = e - bprime
  end function

end module
