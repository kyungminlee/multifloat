module multifloat
  use iso_c_binding

  public :: float64x2


!   type :: t
! contains
!     generic :: operator(+) => t_plus_t, t_plus_i, i_plus_t
!     procedure :: t_plus_t, t_plus_i
!     procedure, pass(rhs) :: i_plus_t
! end type

  type :: float64x2
    real*8 :: limbs(2)
  contains
    generic :: operator(*) => mul, mulf, fmul
    procedure :: mul, mulf
    procedure, pass(b) :: fmul
    generic :: operator(/) => div, divf, fdiv
    procedure :: div, divf
    procedure, pass(b) :: fdiv
  end type

  interface operator (+)
    module procedure add
  end interface

  interface operator (-)
    module procedure sub
  end interface

  interface operator (.lt.)
    module procedure lt
  end interface

  interface assignment(=)
    module procedure assign_from_double
  end interface

  interface operator (.eq.)
    module procedure eq_ff, eq_fd, eq_df
  end interface

  interface operator (.ne.)
    module procedure ne_ff, ne_fd, ne_df
  end interface
contains
  elemental subroutine assign_from_double(lhs, rhs)
    type(float64x2), intent(out) :: lhs
    double precision, intent(in) :: rhs
    lhs%limbs(1) = rhs
    lhs%limbs(2) = 0.0d0
  end subroutine

  pure function eq_ff(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    type(float64x2) :: x2, y2
    x2 = x
    y2 = y
    call renormalize(x2)
    call renormalize(y2)
    z = (x2%limbs(1) .eq. y2%limbs(1)) .and. (x2%limbs(2) .eq. y2%limbs(2))
  end function

  pure function ne_ff(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    z = .not. eq_ff(x, y)
  end function

  pure function eq_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    type(float64x2) :: x2
    x2 = x
    call renormalize(x2)
    z = (x2%limbs(1) .eq. d) .and. (x2%limbs(2) .eq. 0.0d0)
  end function

  pure function eq_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = eq_fd(x, d)
  end function

  pure function ne_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    z = .not. eq_fd(x, d)
  end function

  pure function ne_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = .not. eq_fd(x, d)
  end function

  elemental function div(a, b) result(c)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: c, r
    double precision :: q1, q2
    q1 = a%limbs(1) / b%limbs(1)
    r%limbs(1) = a%limbs(1) - q1 * b%limbs(1)
    r%limbs(2) = a%limbs(2) - q1 * b%limbs(2)
    call fast_two_sum(r%limbs(1), r%limbs(2))
    q2 = r%limbs(1) / b%limbs(1)
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function fdiv(a, b) result(c)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: c
    double precision :: q1, q2, r
    q1 = a / b%limbs(1)
    r = (a - q1 * b%limbs(1)) - q1 * b%limbs(2)
    q2 = r / b%limbs(1)
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function divf(a, b) result(c)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: c
    double precision :: q1, q2, r
    q1 = a%limbs(1) / b
    r = (a%limbs(1) - q1 * b) + a%limbs(2)
    q2 = r / b
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function mul(a, b) result(c)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: c
    c%limbs = a%limbs * b%limbs
  end function

  elemental function fmul(a, b) result(c)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: c
    c%limbs = a * b%limbs
  end function

  elemental function mulf(a, b) result(c)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: c
    c%limbs = a%limbs * b
  end function

  pure subroutine fast_two_sum(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, b_prime, b_err
    s = a + b
    b_prime = s - a
    b_err = b - b_prime
    a = s
    b = b_err
  end

  pure subroutine two_sum(a, b)
    double precision, intent(inout):: a, b
  end subroutine

  elemental subroutine renormalize(x)
    type(float64x2), intent(inout) :: x
    call fast_two_sum(x%limbs(1), x%limbs(2))
  end subroutine

  function add(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: z
    print *, x%limbs
    print *, y%limbs
    z%limbs = x%limbs + y%limbs
    call renormalize(z)
  end function

  elemental function sub(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: z
    z%limbs = x%limbs - y%limbs
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  pure function lt(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: x2, y2
    logical :: z
    x2 = x
    y2 = y
    call renormalize(x2)
    call renormalize(y2)
    z = (x2%limbs(1) .lt. y2%limbs(1)) .or.     &
        & ((x2%limbs(1) .eq. y2%limbs(1)) .and. &
        &  (x2%limbs(2) .lt. y2%limbs(2)))
  end function


end module
