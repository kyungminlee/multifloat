! bind(c) wrappers over the NATIVE float64x2 elemental operators. Used
! by the cross-language differential fuzz driver (test/crosscheck.cc) to
! exercise both the C++ and Fortran implementations of the DD kernels
! that are INDEPENDENTLY implemented on each side — add, sub, mul, div,
! neg, abs, and the ordered/equality comparisons. A bit-for-bit
! divergence on both limbs flags real algorithmic drift between
! src/multifloats_math.cc and fsrc/multifloats.fypp.
!
! Distinct from dd_bindc.f90 in intent: dd_bindc reimplements the DD
! algorithm in pure Fortran to isolate ABI overhead in bench_abi, while
! fnat_* delegates straight to `a + b` / sqrt(a) / etc. so the
! differential compares two production code paths rather than a
! hand-mirrored reimplementation.
module crosscheck_bindings
  ! Import the whole module so the overloaded arithmetic / comparison
  ! operators and the `sqrt` / `abs` generics are visible. An `only:`
  ! list would drop the operator interfaces and the compiler would
  ! reject `type(float64x2) + type(float64x2)`.
  use multifloats
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none
  private

  ! Layout matches float64x2_t in include/multifloats.h and dd_c in
  ! dd_bindc.f90. Declared locally so this module does not depend on
  ! dd_bindc (separate concern).
  type, bind(c), public :: dd_c
    real(c_double) :: hi, lo
  end type dd_c

  public :: fnat_add, fnat_sub, fnat_mul, fnat_div
  public :: fnat_neg, fnat_abs, fnat_sqrt
  public :: fnat_eq,  fnat_ne,  fnat_lt
  ! fnat_le / fnat_ge intentionally omitted: the C++ `<=` is defined as
  ! `!(r < *this)` which is true on NaN pairs (IEEE `<` is false, negation
  ! flips it), whereas Fortran `<=` is IEEE-strict (false on NaN). Any
  ! cross-check would fire on every NaN-bearing iteration for a
  ! pre-existing semantic difference that is out of scope for PR-1.

contains

  pure function c_to_f(c) result(r)
    type(dd_c), intent(in) :: c
    type(float64x2) :: r
    r%limbs(1) = c%hi
    r%limbs(2) = c%lo
  end function

  pure function f_to_c(r) result(c)
    type(float64x2), intent(in) :: r
    type(dd_c) :: c
    c%hi = r%limbs(1)
    c%lo = r%limbs(2)
  end function

  pure function fnat_add(a, b) result(res) bind(c, name='fnat_add')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) + c_to_f(b))
  end function

  pure function fnat_sub(a, b) result(res) bind(c, name='fnat_sub')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) - c_to_f(b))
  end function

  pure function fnat_mul(a, b) result(res) bind(c, name='fnat_mul')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) * c_to_f(b))
  end function

  pure function fnat_div(a, b) result(res) bind(c, name='fnat_div')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) / c_to_f(b))
  end function

  pure function fnat_neg(a) result(res) bind(c, name='fnat_neg')
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    res = f_to_c(-c_to_f(a))
  end function

  pure function fnat_abs(a) result(res) bind(c, name='fnat_abs')
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    res = f_to_c(abs(c_to_f(a)))
  end function

  ! sqrt is included as an ABI smoke test: on the Fortran side sqrt is
  ! delegated to the C++ sqrtdd kernel (see C_DELEGATE_UNARY_MAP in
  ! fsrc/multifloats.fypp), so any mismatch vs. C++ mf::sqrt indicates
  ! ABI corruption rather than algorithmic drift.
  pure function fnat_sqrt(a) result(res) bind(c, name='fnat_sqrt')
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    res = f_to_c(sqrt(c_to_f(a)))
  end function

  pure function fnat_eq(a, b) result(res) bind(c, name='fnat_eq')
    type(dd_c), intent(in), value :: a, b
    integer(c_int) :: res
    res = merge(1, 0, c_to_f(a) == c_to_f(b))
  end function

  pure function fnat_ne(a, b) result(res) bind(c, name='fnat_ne')
    type(dd_c), intent(in), value :: a, b
    integer(c_int) :: res
    res = merge(1, 0, c_to_f(a) /= c_to_f(b))
  end function

  pure function fnat_lt(a, b) result(res) bind(c, name='fnat_lt')
    type(dd_c), intent(in), value :: a, b
    integer(c_int) :: res
    res = merge(1, 0, c_to_f(a) < c_to_f(b))
  end function

end module crosscheck_bindings
