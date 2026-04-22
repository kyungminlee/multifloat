// Complex overloads for the Kind-D parity gap functions — operations the
// C ABI exposes via `c*dd` kernels but which have no std::complex free
// function to specialize (expm1, log2, log10, log1p) or whose name is
// outside the standard library (sinpi, cospi).
//
// Bodies mirror the original extern "C" c*dd implementations in
// multifloats_math_abi_complex.inc; the C shim wrappers over there now
// call these directly.

namespace _abi {
// Local marshaling helpers. The public `to_f64x2` / `from_f64x2` in
// multifloats::detail handle the scalar conversion; wrap the complex
// pair inline since std::complex<float64x2> ↔ complex64x2_t isn't
// worth a separate header declaration.
inline float64x2_t _c2(float64x2 const &x) { return {x._limbs[0], x._limbs[1]}; }
inline float64x2 _f2(float64x2_t x) { float64x2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
inline complex64x2_t _to(std::complex<float64x2> const &z) { return {_c2(z.real()), _c2(z.imag())}; }
inline std::complex<float64x2> _from(complex64x2_t z) { return {_f2(z.re), _f2(z.im)}; }
} // namespace _abi

// log10(z) = log(z) · (1/ln 10). Routes through ::clogdd (the authoritative
// DD complex-log kernel) so the DD precision profile matches the scalar path.
std::complex<float64x2> log10(std::complex<float64x2> const &z) {
  static const float64x2 inv_ln10 = float64x2(0x1.bcb7b1526e50ep-2,
                                              0x1.95355baaafad3p-57);
  std::complex<float64x2> l = _abi::_from(::clogdd(_abi::_to(z)));
  return {l.real() * inv_ln10, l.imag() * inv_ln10};
}

// log2(z) = log(z) · (1/ln 2).
std::complex<float64x2> log2(std::complex<float64x2> const &z) {
  static const float64x2 inv_ln2 =
      float64x2(1.4426950408889634, 2.0355273740931033e-17);
  std::complex<float64x2> l = _abi::_from(::clogdd(_abi::_to(z)));
  return {l.real() * inv_ln2, l.imag() * inv_ln2};
}

// log1p(z) = log(1 + z) without losing precision for |z| << 1. Uses
// log1p on the real scalar form |1+z|² - 1 = 2a + a² + b², which is
// computed without the 1 ± 1 cancellation. The imaginary part is
// atan2(b, 1 + a) — atan2 handles the branch cut and signed-zero
// propagation directly. Matches libquadmath clog1pq.
std::complex<float64x2> log1p(std::complex<float64x2> const &z) {
  float64x2 a = z.real(), b = z.imag();
  float64x2 u = (a + a) + a * a + b * b;          // 2a + a² + b²
  float64x2 half(0.5, 0.0);
  float64x2 r = half * log1p(u);
  float64x2 phi = atan2(b, float64x2(1.0, 0.0) + a);
  return {r, phi};
}

// expm1(z) = e^z - 1 without losing precision for |z| << 1. See the long
// cexpm1dd comment in multifloats_math_abi_complex.inc for the full
// cancellation-guard derivation; this is the same body.
std::complex<float64x2> expm1(std::complex<float64x2> const &z) {
  float64x2 a = z.real(), b = z.imag();
  float64x2 ea_m1 = expm1(a);
  float64x2 ea = exp(a);
  float64x2 s, c;
  sincos(b, s, c);
  float64x2 shalf, chalf;
  sincos(float64x2(0.5, 0.0) * b, shalf, chalf);
  float64x2 cos_m1 = -((shalf * shalf) + (shalf * shalf));
  float64x2 term1 = ea_m1 * c;
  float64x2 re_dd = term1 + cos_m1;
  float64x2 im = ea * s;

  double abs_re    = std::fabs(re_dd._limbs[0]);
  double abs_term1 = std::fabs(term1._limbs[0]);
  double abs_cosm1 = std::fabs(cos_m1._limbs[0]);
  double abs_max   = abs_term1 > abs_cosm1 ? abs_term1 : abs_cosm1;
  constexpr double kCancelThresh = 0x1p-1;
  float64x2 re;
  if (abs_re >= kCancelThresh * abs_max || abs_max == 0.0) {
    re = re_dd;
  } else {
    detail::float64x3 ea_td = exp_full_td(a);
    detail::float64x3 s_td, c_td;
    sincos_td(b, s_td, c_td);
    detail::float64x3 prod_td = detail::td_mul_td(ea_td, c_td);
    detail::float64x3 diff_td = detail::td_sub_double(prod_td, 1.0);
    re = detail::td_to_dd(diff_td);
  }
  return {re, im};
}

// sin(π·z) = sin(π·a) cosh(π·b) + i cos(π·a) sinh(π·b).
std::complex<float64x2> sinpi(std::complex<float64x2> const &z) {
  static const float64x2 pi_dd(pi_dd_hi, pi_dd_lo);
  float64x2 a = z.real(), b = z.imag();
  float64x2 s = sinpi(a);
  float64x2 c = cospi(a);
  float64x2 pib = pi_dd * b;
  float64x2 sh, ch;
  sinhcosh(pib, sh, ch);
  return {s * ch, c * sh};
}

// cos(π·z) = cos(π·a) cosh(π·b) − i sin(π·a) sinh(π·b).
std::complex<float64x2> cospi(std::complex<float64x2> const &z) {
  static const float64x2 pi_dd(pi_dd_hi, pi_dd_lo);
  float64x2 a = z.real(), b = z.imag();
  float64x2 s = sinpi(a);
  float64x2 c = cospi(a);
  float64x2 pib = pi_dd * b;
  float64x2 sh, ch;
  sinhcosh(pib, sh, ch);
  return {c * ch, -(s * sh)};
}
