module multifloat
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic

  public :: float64x2
  public :: abs, sqrt, sign, min, max, floor, ceiling, aint, anint, exp, log, log10

  type :: float64x2
    real*8 :: limbs(2)
  contains
    generic :: operator(*) => mul, mulf, fmul
    procedure :: mul, mulf
    procedure, pass(b) :: fmul
    generic :: operator(/) => div, divf, fdiv
    procedure :: div, divf
    procedure, pass(b) :: fdiv
    procedure :: isfinite => isfinite_f
    procedure :: precision => precision_f
    procedure :: minexponent => minexponent_f
    procedure :: maxexponent => maxexponent_f
    procedure :: tiny => tiny_f
    procedure :: huge => huge_f
    procedure :: exponent => exponent_f
    procedure :: scale => scale_f
    procedure :: floor => floor_f
    procedure :: ceiling => ceiling_f
    procedure :: aint => aint_f
    procedure :: anint => anint_f
    procedure :: abs => abs_f
    procedure :: sqrt => sqrt_f
    procedure :: fraction => fraction_f
    procedure :: set_exponent => set_exponent_f
    procedure :: spacing => spacing_f
    procedure :: rrspacing => rrspacing_f
    procedure :: exp => exp_f
    procedure :: log => log_f
    procedure :: log10 => log10_f
  end type

  interface operator (+)
    module procedure add, add_fd, add_df
  end interface

  interface operator (-)
    module procedure sub, sub_fd, sub_df, neg
  end interface

  interface operator (.lt.)
    module procedure lt, lt_fd, lt_df
  end interface

  interface operator (.gt.)
    module procedure gt, gt_fd, gt_df
  end interface

  interface operator (.le.)
    module procedure le, le_fd, le_df
  end interface

  interface operator (.ge.)
    module procedure ge, ge_fd, ge_df
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

  interface abs
    module procedure abs_f
  end interface

  interface sqrt
    module procedure sqrt_f
  end interface

  interface sign
    module procedure sign_f, sign_fd, sign_df
  end interface

  interface min
    module procedure min_ff, min_fd, min_df
  end interface

  interface max
    module procedure max_ff, max_fd, max_df
  end interface

  interface aint
    module procedure aint_f
  end interface

  interface anint
    module procedure anint_f
  end interface

  interface floor
    module procedure floor_f
  end interface

  interface ceiling
    module procedure ceiling_f
  end interface

  interface exp
    module procedure exp_f
  end interface

  interface exp2
    module procedure exp2_f
  end interface

  interface log
    module procedure log_f
  end interface

  interface log10
    module procedure log10_f
  end interface

  ! Constants for bounds
  real(dp), parameter :: EXP2_MIN = -1022.0_dp ! approx -0x1.FF00000000000p+0009
  real(dp), parameter :: EXP2_MAX =  1024.0_dp ! approx +0x1.FFFFFFFFFFFFFp+0009

  ! Log2(e) as a double-double
  ! hi = +0x1.71547652B82FEp+0000, lo = +0x1.777D0FFDA0D24p-0056
  type(float64x2), parameter :: LOG2_E = float64x2((/dble(z'3ff71547652b82fe'), dble(z'3c7777d0ffda0d24')/))

  ! LN_2 = (+0x1.62E42FEFA39EFp-0001, +0x1.ABC9E3B39803Fp-0056)
  type(float64x2), parameter :: LN_2 = float64x2([ &
    0.69314718055994529_dp,  &
    2.3190468138462996e-17_dp ])

  ! LOG10_2 = (+0x1.34413509F79FFp-0002, -0x1.9DC1DA994FD21p-0059)
  type(float64x2), parameter :: LOG10_2 = float64x2([ &
    0.30102999566398119_dp, &
    -2.8037281277851703e-18_dp ])

contains

  function exp_f(x) result(res)
    ! class(float64x2), intent(in) :: x
    ! type(float64x2) :: z
    ! real(16) :: q
    ! if (.not. ieee_is_finite(x%limbs(1))) then
    !    z%limbs(1) = exp(x%limbs(1))
    !    z%limbs(2) = 0.0d0
    !    return
    ! end if
    ! q = exp(real(x%limbs(1), 16) + real(x%limbs(2), 16))
    ! z%limbs(1) = dble(q)
    ! z%limbs(2) = dble(q - real(z%limbs(1), 16))
    implicit none
    class(float64x2), intent(in) :: x
    type(float64x2)             :: res
    type(float64x2)             :: x_log2e

    ! Requires a double-double multiplication function (mul)
    x_log2e = mul(x, LOG2_E)
    res = exp2_f(x_log2e)
  end function

  function exp2_f(x) result(res)
    use, intrinsic :: ieee_arithmetic
    implicit none
    class(float64x2), intent(in) :: x
    type(float64x2)             :: res
    type(float64x2)             :: r, p
    real(dp)                  :: n_float
    integer                   :: n_int, half_n
    ! print *, "Checking ", x%limbs(1)

    ! 1. Check Overflow / Underflow limits (based on the hi limb)
    if (x%limbs(1) < EXP2_MIN) then
      ! print *, x%limbs(1), " Less than min"
      res = float64x2((/0.0_dp, 0.0_dp/))
      return
    else if (x%limbs(1) > EXP2_MAX) then
      ! print *, x%limbs(1), " Greater than max"
      ! Return Infinity (using standard IEEE positive infinity generation)
      res = float64x2((/ieee_value(1.0d0, ieee_positive_inf), 0.0_dp/))
      ! res = float64x2((/huge(1.0_dp) * 2.0_dp, 0.0_dp/))
      return
    else
      ! print *, x%limbs(1), " Is OK"
    end if

    ! 2. Range Reduction: Separate into integer and fractional parts
    ! anint() is equivalent to Julia's trunc(first(x) + copysign(0.5, first(x)))
    n_float = anint(x%limbs(1)) 

    ! r = x - n_float (Requires double-double addition)
    r = add(x, float64x2((/-n_float, 0.0_dp/)))

    ! Scale r by 1/8 (0.125)
    r%limbs(1) = r%limbs(1) * 0.125_dp
    r%limbs(2) = r%limbs(2) * 0.125_dp

    ! print *, "A: ", r%limbs
    !
    ! 3. Polynomial Evaluation
    p = evaluate_exp2_poly(r)

    ! print *, "p: ", p
    !
    ! 4. Undo the 1/8 scaling by squaring the result 3 times (p^8)
    p = mul(p,p)
    p = mul(p,p)
    p = mul(p,p)

    ! print *, "p2: ", p
    !
    ! 5. Scale by 2^n (Undo the integer part reduction)
    ! Split n into two halves to prevent intermediate overflow before the end
    n_int = int(n_float)
    half_n = shifta(n_int, 1) ! equivalent to n >> 1

    ! Fortran's intrinsic `scale(x, i)` is exactly C's `ldexp`, multiplying x by 2^i
    p%limbs(1) = scale(p%limbs(1), half_n)
    p%limbs(2) = scale(p%limbs(2), half_n)

    ! print *, "half_n: ", half_n
    ! print *, "p3: ", p
    !
    p%limbs(1) = scale(p%limbs(1), n_int - half_n)
    p%limbs(2) = scale(p%limbs(2), n_int - half_n)

    ! print *, "p4: ", p
    !
    !
    res = p
  end function exp2_f

  elemental function evaluate_exp2_poly(x) result(res)
    class(float64x2), intent(in) :: x
    type(float64x2)             :: res

    ! These constants are exactly mapped from Julia's `_exp2_coefficients(Float64, Val{2})`
    ! Note: Converting hex-floats to decimal constants ensures cross-compiler Fortran safety.

    ! ! c13 = (+0x1.816519F74C4AFp-0040)
    ! res = float64x2((/1.3683838383401736e-12_dp, 0.0_dp/))
    !
    ! ! c12 = (+0x1.C3C1919538484p-0036)
    ! res = add(mul(res, x), float64x2((/2.5694246028169222e-11_dp, 0.0_dp/)))
    !
    ! ! c11 = (+0x1.E8CAC72F6E9E5p-0032)
    ! res = add(mul(res, x), float64x2((/4.4449852509176182e-10_dp, 0.0_dp/)))
    !
    ! ! c10 = (+0x1.E4CF5152FBB30p-0028)
    ! res = add(mul(res, x), float64x2((/7.0700589886981883e-09_dp, 0.0_dp/)))
    !
    ! ! c9  = (+0x1.B5253D395E80Fp-0024)
    ! res = add(mul(res, x), float64x2((/1.0182436893699479e-07_dp, 0.0_dp/)))
    !
    ! ! c8  = (+0x1.62C0223A5C863p-0020, -0x1.99EF542AA8E1Ep-0074)
    ! res = add(mul(res, x), float64x2((/1.3215486790144309e-06_dp, -8.4552438848148008e-23_dp/)))
    !
    ! ! c7  = (+0x1.FFCBFC588B0C7p-0017, -0x1.E645E286FE571p-0071)
    ! res = add(mul(res, x), float64x2((/1.5252733804059837e-05_dp, -8.0558137352329618e-22_dp/)))
    !
    ! ! c6  = (+0x1.430912F86C787p-0013, +0x1.BC7CDBCDC0339p-0067)
    ! res = add(mul(res, x), float64x2((/0.00015403530393381609_dp, 1.1764353493649931e-20_dp/)))
    !
    ! ! c5  = (+0x1.5D87FE78A6731p-0010, +0x1.0717F88815ADFp-0066)
    ! res = add(mul(res, x), float64x2((/0.0013333558146428443_dp,  1.3934305888941271e-20_dp/)))
    !
    ! ! c4  = (+0x1.3B2AB6FBA4E77p-0007, +0x1.4E65DFEF67D34p-0062)
    ! res = add(mul(res, x), float64x2((/0.0096181291076284762_dp,  2.8368536120853503e-19_dp/)))
    !
    ! ! c3  = (+0x1.C6B08D704A0C0p-0005, -0x1.D3316275139AEp-0059)
    ! res = add(mul(res, x), float64x2((/0.055504108664821578_dp,  -3.1793575936306567e-18_dp/)))
    !
    ! ! c2  = (+0x1.EBFBDFF82C58Fp-0003, -0x1.5E43A53E454F1p-0057)
    ! res = add(mul(res, x), float64x2((/0.24022650695910069_dp,   -9.5080088927063464e-18_dp/)))
    !
    ! ! c1  = (+0x1.62E42FEFA39EFp-0001, +0x1.ABC9E3B39803Fp-0056)
    ! res = add(mul(res, x), float64x2((/0.69314718055994529_dp,   2.3190468138462996e-17_dp/)))
    !
    ! ! c0  = (+0x1.0000000000000p+0000, +0x1.314BACF0323FFp-0113)
    ! res = add(mul(res, x), float64x2((/1.0_dp,                   1.1491122822830842e-34_dp/)))


    ! 3d7816519f74c4af      0
    ! 4429314773342274735   0
    ! 1.36919783052689e-12  0.0
    res = float64x2((/dble(z'3d7816519f74c4af'), dble(z'0')/))
    ! i=13
    ! 3dbc3c1919538484      0
    ! 4448496610431960196   0
    ! 2.567936278800348e-11 0.0
    res = add(mul(res, x), float64x2((/dble(z'3dbc3c1919538484'), dble(z'0')/)))



    ! 3dfe8cac72f6e9e5      0
    res = add(mul(res, x), float64x2((/dble(z'3dfe8cac72f6e9e5'), dble(z'0')/)))
    ! 3e3e4cf5152fbb30      0
    res = add(mul(res, x), float64x2((/dble(z'3e3e4cf5152fbb30'), dble(z'0')/)))
    ! 3e7b5253d395e80f      0
    res = add(mul(res, x), float64x2((/dble(z'3e7b5253d395e80f'), dble(z'0')/)))
    ! 3eb62c0223a5c863      bb599ef542aa8e1e
    res = add(mul(res, x), float64x2((/dble(z'3eb62c0223a5c863'), dble(z'bb599ef542aa8e1e')/)))
    ! 3eeffcbfc588b0c7      bb8e645e286fe571
    res = add(mul(res, x), float64x2((/dble(z'3eeffcbfc588b0c7'), dble(z'bb8e645e286fe571')/)))
    ! 3f2430912f86c787      3bcbc7cdbcdc0339
    res = add(mul(res, x), float64x2((/dble(z'3f2430912f86c787'), dble(z'3bcbc7cdbcdc0339')/)))
    ! 3f55d87fe78a6731      3bd0717f88815adf
    res = add(mul(res, x), float64x2((/dble(z'3f55d87fe78a6731'), dble(z'3bd0717f88815adf')/)))
    ! 3f83b2ab6fba4e77      3c14e65dfef67d34
    res = add(mul(res, x), float64x2((/dble(z'3f83b2ab6fba4e77'), dble(z'3c14e65dfef67d34')/)))
    ! 3fac6b08d704a0c0      bc4d3316275139ae
    res = add(mul(res, x), float64x2((/dble(z'3fac6b08d704a0c0'), dble(z'bc4d3316275139ae')/)))
    ! 3fcebfbdff82c58f      bc65e43a53e454f1
    res = add(mul(res, x), float64x2((/dble(z'3fcebfbdff82c58f'), dble(z'bc65e43a53e454f1')/)))
    ! 3fe62e42fefa39ef      3c7abc9e3b39803f
    res = add(mul(res, x), float64x2((/dble(z'3fe62e42fefa39ef'), dble(z'3c7abc9e3b39803f')/)))
    ! 3ff0000000000000      38e314bacf0323ff
    res = add(mul(res, x), float64x2((/dble(z'3ff0000000000000'), dble(z'38e314bacf0323ff')/)))
  end function evaluate_exp2_poly


  ! ! =========================================================================
  ! ! Base-e Logarithm: log(x)
  ! ! =========================================================================
  ! function f64x2_log(x) result(res)
  !   type(float64x2), intent(in) :: x
  !   type(float64x2)             :: res
  !
  !   ! 1. Check special limits
  !   if (x%limbs(1) < 0.0_dp .or. x%limbs(1) /= x%limbs(1)) then
  !     res = float64x2([ieee_value(1.0_dp, ieee_quiet_nan), 0.0_dp])
  !     return
  !   else if (x%limbs(1) == 0.0_dp) then
  !     res = float64x2([-huge(1.0_dp) * 2.0_dp, 0.0_dp]) ! -Infinity
  !     return
  !   else if (x%limbs(1) > huge(1.0_dp)) then
  !     res = x ! +Infinity
  !     return
  !   end if
  !
  !   ! 2. Compute log2(x) * ln(2)
  !   res = mul(f64x2_log2(x), LN_2)
  ! end function f64x2_log
  !
  !
  ! ! =========================================================================
  ! ! Base-10 Logarithm: log10(x)
  ! ! =========================================================================
  ! function f64x2_log10(x) result(res)
  !   type(float64x2), intent(in) :: x
  !   type(float64x2)             :: res
  !
  !   ! 1. Check special limits
  !   if (x%limbs(1) < 0.0_dp .or. x%limbs(1) /= x%limbs(1)) then
  !     res = float64x2([ieee_value(1.0_dp, ieee_quiet_nan), 0.0_dp])
  !     return
  !   else if (x%limbs(1) == 0.0_dp) then
  !     res = float64x2([-huge(1.0_dp) * 2.0_dp, 0.0_dp]) ! -Infinity
  !     return
  !   else if (x%limbs(1) > huge(1.0_dp)) then
  !     res = x ! +Infinity
  !     return
  !   end if
  !
  !   ! 2. Compute log2(x) * log10(2)
  !   res = mul(f64x2_log2(x), LOG10_2)
  ! end function f64x2_log10
  !
  ! ! =========================================================================
  ! ! Core Base-2 Logarithm: unsafe_log2(x)
  ! ! =========================================================================
  ! function f64x2_log2(x) result(res)
  !   type(float64x2), intent(in) :: x
  !   type(float64x2)             :: res
  !   type(float64x2)             :: m, t_direct, p_direct, t_table, p_table, t_sqr
  !   type(float64x2)             :: center, value
  !   real(dp)                    :: e_float, direct_lo, direct_hi
  !   integer                     :: e_int, index
  !
  !   ! Limits for the direct polynomial evaluation (no table lookup needed)
  !   direct_lo = 15.0_dp / 16.0_dp
  !   direct_hi = 17.0_dp / 16.0_dp
  !
  !   ! Standard identity for x
  !   type(float64x2), parameter :: ONE = float64x2([1.0_dp, 0.0_dp])
  !
  !   if (x%limbs(1) > direct_lo .and. x%limbs(1) < direct_hi) then
  !     ! --- Fast Path: Near 1.0, use the wide polynomial ---
  !     ! t_direct = (x - 1) / (x + 1)
  !     t_direct = div(sub(x, ONE), add(x, ONE))
  !     t_sqr    = mul(t_direct, t_direct)
  !     p_direct = evaluate_log2_wide(t_sqr)
  !
  !     res = mul(t_direct, p_direct)
  !     return
  !   else
  !     ! --- Table-Assisted Path ---
  !     ! Extract exponent and fractional part. `fraction` instrinsic is 
  !     ! equivalent to `m = unsafe_ldexp(x, -e)`
  !     ! Note: Fortran's `fraction(x)` maps strictly to [0.5, 1.0)
  !     e_int   = exponent(x%limbs(1))
  !     e_float = real(e_int, dp)
  !
  !     ! Scale down x to get the mantissa
  !     m%limbs(1) = scale(x%limbs(1), -e_int)
  !     m%limbs(2) = scale(x%limbs(2), -e_int)
  !
  !     ! Get table index using top 5 bits of mantissa
  !     ! (Implementation of `_log2_table_index` omitted for brevity; it requires
  !     !  extracting bits 47-51 of the IEEE-754 mantissa to fetch index 1-32).
  !     index = get_table_index(m%limbs(1))
  !
  !     ! Fetch precomputed log2(center) lookup
  !     call get_log2_table_values(index, center, value)
  !
  !     ! t_table = (m - center) / (m + center)
  !     t_table = div(sub(m, center), add(m, center))
  !     t_sqr   = mul(t_table, t_table)
  !     p_table = evaluate_log2_narrow(t_sqr)
  !
  !     ! result = e + value + (t_table * p_table)
  !     res = add(float64x2([e_float, 0.0_dp]), value)
  !     res = add(res, mul(t_table, p_table))
  !   end if
  ! end function f64x2_log2
  !
  ! ! =========================================================================
  ! ! Polynomial Evaluators (Horner's Method)
  ! ! =========================================================================
  ! function evaluate_log2_wide(x) result(res)
  !   type(float64x2), intent(in) :: x
  !   type(float64x2)             :: res
  !
  !   ! Maps to _log2_kernel_coefficients_wide(Float64, Val{2})
  !   ! c8 = (+0x1.5D108CD7E21EBp-0003,)
  !   res = float64x2([ 0.17044199464539825_dp, 0.0_dp ])
  !
  !   ! c7 = (+0x1.89F2F69137330p-0003,)
  !   res = add(mul(res, x), float64x2([ 0.19235779058694038_dp, 0.0_dp ]))
  !
  !   ! ... remaining coefficients up to c0 inserted here using `add(mul(res, x), c_n)`...
  !
  !   ! c0 = (+0x1.71547652B82FEp+0001, +0x1.777D0FFDA0D24p-0055)
  !   res = add(mul(res, x), float64x2([ 2.8853900817779268_dp, 4.0863032542385496e-17_dp ]))
  ! end function evaluate_log2_wide
  !
  ! function evaluate_log2_narrow(x) result(res)
  !   type(float64x2), intent(in) :: x
  !   type(float64x2)             :: res
  !
  !   ! Maps to _log2_kernel_coefficients_narrow(Float64, Val{2})
  !   ! c6 = (+0x1.C6A48D52BA6C7p-0003,)
  !   res = float64x2([ 0.22199347895287959_dp, 0.0_dp ])
  !
  !   ! ... remaining coefficients up to c0 inserted here ...
  !
  !   ! c0 = (+0x1.71547652B82FEp+0001, +0x1.777D0FFDA0D24p-0055)
  !   res = add(mul(res, x), float64x2([ 2.8853900817779268_dp, 4.0863032542385496e-17_dp ]))
  ! end function evaluate_log2_narrow
  !
  ! ! (Stubs for get_table_index and get_log2_table_values go here)
  !
  !
  !



  elemental function log_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    real(16) :: q
    if (.not. ieee_is_finite(x%limbs(1)) .or. x%limbs(1) <= 0.0d0) then
       z%limbs(1) = log(x%limbs(1))
       z%limbs(2) = 0.0d0
       return
    end if
    q = log(real(x%limbs(1), 16) + real(x%limbs(2), 16))
    z%limbs(1) = dble(q)
    z%limbs(2) = dble(q - real(z%limbs(1), 16))
  end function

  elemental function log10_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    real(16) :: q
    if (.not. ieee_is_finite(x%limbs(1)) .or. x%limbs(1) <= 0.0d0) then
       z%limbs(1) = log10(x%limbs(1))
       z%limbs(2) = 0.0d0
       return
    end if
    q = log10(real(x%limbs(1), 16) + real(x%limbs(2), 16))
    z%limbs(1) = dble(q)
    z%limbs(2) = dble(q - real(z%limbs(1), 16))
  end function

  elemental function to_f64x2_d(d) result(z)
    double precision, intent(in) :: d
    type(float64x2) :: z
    z%limbs(1) = d
    z%limbs(2) = 0.0d0
  end function

  elemental function abs_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    if (x%limbs(1) < 0.0d0) then
      z = -x
    else
      z = x
    end if
  end function

  elemental function sqrt_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    double precision :: s, p, e, h
    if (x%limbs(1) <= 0.0d0) then
       if (x%limbs(1) < 0.0d0) then
          z%limbs(1) = ieee_value(0.0d0, ieee_quiet_nan)
          z%limbs(2) = 0.0d0
       else
          z%limbs = 0.0d0
       end if
       return
    end if
    s = sqrt(x%limbs(1))
    ! Refinement: z = s + (x - s*s) / (2*s)
    call two_prod(s, s, p, e)
    h = (x%limbs(1) - p) + (x%limbs(2) - e)
    ! Must use double-double addition
    z%limbs(1) = s
    z%limbs(2) = h / (2.0d0 * s)
    call renormalize(z)
  end function

  elemental function sign_f(a, b) result(z)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: z
    z = abs(a)
    if (b%limbs(1) < 0.0d0) z = -z
  end function

  elemental function sign_fd(a, b) result(z)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: z
    z = abs(a)
    if (b < 0.0d0) z = -z
  end function

  elemental function sign_df(a, b) result(z)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: z
    z%limbs(1) = abs(a)
    z%limbs(2) = 0.0d0
    if (b%limbs(1) < 0.0d0) z = -z
  end function

  elemental function min_ff(a, b) result(z)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: z
    if (a < b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function min_fd(a, b) result(z)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: z
    if (a < b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function min_df(a, b) result(z)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: z
    if (a < b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function max_ff(a, b) result(z)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: z
    if (a > b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function max_fd(a, b) result(z)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: z
    if (a > b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function max_df(a, b) result(z)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: z
    if (a > b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function aint_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = aint(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = aint(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function anint_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = anint(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = anint(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function floor_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = floor(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = floor(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function ceiling_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = ceiling(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = ceiling(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function precision_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = 31
  end function

  elemental function minexponent_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = minexponent(1.0d0)
  end function

  elemental function maxexponent_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = maxexponent(1.0d0)
  end function

  elemental function tiny_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = tiny(1.0d0)
    z%limbs(2) = 0.0d0
  end function

  elemental function huge_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = huge(1.0d0)
    z%limbs(2) = 0.0d0
  end function

  elemental function exponent_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = exponent(x%limbs(1))
  end function

  elemental function scale_f(x, i) result(z)
    class(float64x2), intent(in) :: x
    integer, intent(in) :: i
    type(float64x2) :: z
    z%limbs = scale(x%limbs, i)
  end function

  elemental function fraction_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = scale_f(x, -exponent(x%limbs(1)))
  end function

  elemental function set_exponent_f(x, i) result(z)
    class(float64x2), intent(in) :: x
    integer, intent(in) :: i
    type(float64x2) :: z
    z = scale_f(fraction_f(x), i)
  end function

  elemental function spacing_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = scale(spacing(x%limbs(1)), -53)
    z%limbs(2) = 0.0d0
  end function

  elemental function rrspacing_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = abs(x) / spacing_f(x)
  end function

  elemental function neg(x) result(z)
    type(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs = -x%limbs
  end function

  elemental function isfinite_f(x) result(z)
    class(float64x2), intent(in) :: x
    logical :: z
    z = ieee_is_finite(x%limbs(1))
  end function

  elemental subroutine assign_from_double(lhs, rhs)
    type(float64x2), intent(out) :: lhs
    double precision, intent(in) :: rhs
    lhs%limbs(1) = rhs
    lhs%limbs(2) = 0.0d0
  end subroutine

  pure function eq_ff(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .eq. y%limbs(1))
      return
    end if
    z = (x%limbs(1) .eq. y%limbs(1)) .and. (x%limbs(2) .eq. y%limbs(2))
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
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .eq. d)
      return
    end if
    z = (x%limbs(1) .eq. d) .and. (x%limbs(2) .eq. 0.0d0)
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
    type(float64x2) :: c
    double precision :: q1, q2, r
    if (b%limbs(1) == 0.0d0) then
      c%limbs(1) = a%limbs(1) / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a%limbs(1) / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    q1 = a%limbs(1) / b%limbs(1)
    r = ieee_fma(-q1, b%limbs(1), a%limbs(1))
    r = r + a%limbs(2) - q1 * b%limbs(2)
    q2 = r / b%limbs(1)
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function fdiv(a, b) result(c)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: c
    double precision :: q1, q2, r
    if (b%limbs(1) == 0.0d0) then
      c%limbs(1) = a / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    q1 = a / b%limbs(1)
    r = ieee_fma(-q1, b%limbs(1), a)
    r = r - q1 * b%limbs(2)
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
    if (b == 0.0d0) then
      c%limbs(1) = a%limbs(1) / b
      c%limbs(2) = 0.0d0
      return
    end if
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b)) then
      c%limbs(1) = a%limbs(1) / b
      c%limbs(2) = 0.0d0
      return
    end if
    q1 = a%limbs(1) / b
    r = ieee_fma(-q1, b, a%limbs(1))
    r = r + a%limbs(2)
    q2 = r / b
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function mul(a, b) result(c)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: c
    double precision :: p, e
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a%limbs(1) * b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    call two_prod(a%limbs(1), b%limbs(1), p, e)
    e = e + a%limbs(1) * b%limbs(2)
    e = e + a%limbs(2) * b%limbs(1)
    e = e + a%limbs(2) * b%limbs(2)
    c%limbs(1) = p
    c%limbs(2) = e
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function fmul(a, b) result(c)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: c
    double precision :: p, e
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a * b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    call two_prod(a, b%limbs(1), p, e)
    e = e + a * b%limbs(2)
    c%limbs(1) = p
    c%limbs(2) = e
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function mulf(a, b) result(c)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: c
    double precision :: p, e
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b)) then
      c%limbs(1) = a%limbs(1) * b
      c%limbs(2) = 0.0d0
      return
    end if
    call two_prod(a%limbs(1), b, p, e)
    e = e + a%limbs(2) * b
    c%limbs(1) = p
    c%limbs(2) = e
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental subroutine fast_two_sum(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, b_prime, b_err
    s = a + b
    b_prime = s - a
    b_err = b - b_prime
    a = s
    b = b_err
  end subroutine

  elemental subroutine two_sum(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, a_prime, b_prime, a_err, b_err
    s = a + b
    a_prime = s - b
    b_prime = s - a_prime
    a_err = a - a_prime
    b_err = b - b_prime
    a = s
    b = a_err + b_err
  end subroutine

  elemental subroutine two_prod(a, b, p, e)
    double precision, intent(in) :: a, b
    double precision, intent(out) :: p, e
    p = a * b
    e = ieee_fma(a, b, -p)
  end subroutine

  elemental subroutine renormalize(x)
    type(float64x2), intent(inout) :: x
    if (.not. ieee_is_finite(x%limbs(1))) return
    call two_sum(x%limbs(1), x%limbs(2))
  end subroutine

  elemental function add(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z%limbs(1) = x%limbs(1) + y%limbs(1)
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = y%limbs(1)
    call two_sum(s, e)
    e = e + x%limbs(2) + y%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function add_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z%limbs(1) = x%limbs(1) + d
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = d
    call two_sum(s, e)
    e = e + x%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function add_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = add_fd(x, d)
  end function

  elemental function sub(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z%limbs(1) = x%limbs(1) - y%limbs(1)
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = -y%limbs(1)
    call two_sum(s, e)
    e = e + x%limbs(2) - y%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function sub_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z%limbs(1) = x%limbs(1) - d
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = -d
    call two_sum(s, e)
    e = e + x%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function sub_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z%limbs(1) = d - x%limbs(1)
      z%limbs(2) = 0.0d0
      return
    end if
    s = d
    e = -x%limbs(1)
    call two_sum(s, e)
    e = e - x%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  pure function lt(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: x2, y2
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .lt. y%limbs(1))
      return
    end if
    x2 = x
    y2 = y
    call renormalize(x2)
    call renormalize(y2)
    z = (x2%limbs(1) .lt. y2%limbs(1)) .or.     &
        & ((x2%limbs(1) .eq. y2%limbs(1)) .and. &
        &  (x2%limbs(2) .lt. y2%limbs(2)))
  end function

  pure function lt_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    type(float64x2) :: x2
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .lt. d)
      return
    end if
    x2 = x
    call renormalize(x2)
    z = (x2%limbs(1) .lt. d) .or. &
        & ((x2%limbs(1) .eq. d) .and. (x2%limbs(2) .lt. 0.0d0))
  end function

  pure function lt_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    type(float64x2) :: x2
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z = (d .lt. x%limbs(1))
      return
    end if
    x2 = x
    call renormalize(x2)
    z = (d .lt. x2%limbs(1)) .or. &
        & ((d .eq. x2%limbs(1)) .and. (0.0d0 .lt. x2%limbs(2)))
  end function

  pure function gt(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    z = lt(y, x)
  end function

  pure function gt_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    z = lt_df(d, x)
  end function

  pure function gt_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = lt_fd(x, d)
  end function

  pure function le(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .le. y%limbs(1))
      return
    end if
    z = .not. lt(y, x)
  end function

  pure function le_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .le. d)
      return
    end if
    z = .not. lt_df(d, x)
  end function

  pure function le_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z = (d .le. x%limbs(1))
      return
    end if
    z = .not. lt_fd(x, d)
  end function

  pure function ge(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .ge. y%limbs(1))
      return
    end if
    z = .not. lt(x, y)
  end function

  pure function ge_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .ge. d)
      return
    end if
    z = .not. lt_fd(x, d)
  end function

  pure function ge_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z = (d .ge. x%limbs(1))
      return
    end if
    z = .not. lt_df(d, x)
  end function

end module
