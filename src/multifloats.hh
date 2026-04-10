#pragma once

#include <cmath>
#include <cstddef>
#include <limits>

namespace multifloats {

namespace detail {

// =============================================================================
// Error-free transformations (inputs by value to avoid aliasing pitfalls)
// =============================================================================

template <typename T>
constexpr void two_sum(T a, T b, T &sum, T &err) {
  T s = a + b;
  T a_prime = s - b;
  T b_prime = s - a_prime;
  err = (a - a_prime) + (b - b_prime);
  sum = s;
}

template <typename T>
constexpr void fast_two_sum(T a, T b, T &sum, T &err) {
  T s = a + b;
  T b_prime = s - a;
  err = b - b_prime;
  sum = s;
}

template <typename T> constexpr T one_prod(T a, T b) { return a * b; }

template <typename T>
constexpr void two_prod(T a, T b, T &prod, T &err) {
  T p = a * b;
  err = std::fma(a, b, -p);
  prod = p;
}

} // namespace detail

template <typename T, std::size_t N> struct MultiFloat {
  static_assert(N == 1 || N == 2, "only N = 1 and N = 2 are implemented");
  T _limbs[N] = {};

  constexpr MultiFloat() = default;
  constexpr MultiFloat(MultiFloat const &) = default;
  constexpr MultiFloat(MultiFloat &&) = default;
  constexpr MultiFloat &operator=(MultiFloat const &) = default;
  constexpr MultiFloat &operator=(MultiFloat &&) = default;

  constexpr MultiFloat(T const &arg) { _limbs[0] = arg; }

  constexpr explicit operator T() const { return _limbs[0]; }

  constexpr bool operator==(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (!(_limbs[i] == rhs._limbs[i])) {
        return false;
      }
    }
    return true;
  }

  constexpr bool operator!=(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] != rhs._limbs[i]) {
        return true;
      }
    }
    return false;
  }

  constexpr bool operator<(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] < rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] < _limbs[i]) {
        return false;
      }
    }
    return false;
  }

  constexpr bool operator>(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] > rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] > _limbs[i]) {
        return false;
      }
    }
    return false;
  }

  constexpr bool operator<=(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] < rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] < _limbs[i]) {
        return false;
      }
    }
    return true;
  }

  constexpr bool operator>=(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] > rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] > _limbs[i]) {
        return false;
      }
    }
    return true;
  }

  constexpr MultiFloat operator+() const { return *this; }

  constexpr MultiFloat operator-() const {
    MultiFloat r;
    for (std::size_t i = 0; i < N; ++i) {
      r._limbs[i] = -_limbs[i];
    }
    return r;
  }

  // ---------------------------------------------------------------------------
  // Binary arithmetic — kernels inlined directly, translated from
  // MultiFloats.jl (mfadd / mfmul) and the Float64x2 division kernel in
  // fsrc/multifloats.fypp.
  // ---------------------------------------------------------------------------

  constexpr MultiFloat operator+(MultiFloat const &rhs) const {
    MultiFloat out;
    if constexpr (N == 1) {
      out._limbs[0] = _limbs[0] + rhs._limbs[0];
    } else { // N == 2
      T a, b, c, d;
      detail::two_sum(_limbs[0], rhs._limbs[0], a, b);
      detail::two_sum(_limbs[1], rhs._limbs[1], c, d);
      detail::fast_two_sum(a, c, a, c);
      b += d;
      b += c;
      detail::fast_two_sum(a, b, out._limbs[0], out._limbs[1]);
    }
    return out;
  }

  constexpr MultiFloat operator-(MultiFloat const &rhs) const {
    return (*this) + (-rhs);
  }

  constexpr MultiFloat operator*(MultiFloat const &rhs) const {
    MultiFloat out;
    if constexpr (N == 1) {
      out._limbs[0] = _limbs[0] * rhs._limbs[0];
    } else { // N == 2
      T p00, e00;
      detail::two_prod(_limbs[0], rhs._limbs[0], p00, e00);
      T p01 = detail::one_prod(_limbs[0], rhs._limbs[1]);
      T p10 = detail::one_prod(_limbs[1], rhs._limbs[0]);
      p01 += p10;
      e00 += p01;
      detail::fast_two_sum(p00, e00, out._limbs[0], out._limbs[1]);
    }
    return out;
  }

  constexpr MultiFloat operator/(MultiFloat const &rhs) const {
    if constexpr (N == 1) {
      MultiFloat out;
      out._limbs[0] = _limbs[0] / rhs._limbs[0];
      return out;
    } else { // N == 2 — single Newton refinement
      MultiFloat u;
      u._limbs[0] = T(1) / rhs._limbs[0];
      MultiFloat quotient = (*this) * u;
      MultiFloat residual = quotient * rhs - (*this);
      MultiFloat correction = residual * u;
      return quotient - correction;
    }
  }

  constexpr MultiFloat &operator+=(MultiFloat const &rhs) {
    return *this = *this + rhs;
  }
  constexpr MultiFloat &operator-=(MultiFloat const &rhs) {
    return *this = *this - rhs;
  }
  constexpr MultiFloat &operator*=(MultiFloat const &rhs) {
    return *this = *this * rhs;
  }
  constexpr MultiFloat &operator/=(MultiFloat const &rhs) {
    return *this = *this / rhs;
  }
};

// =============================================================================
// <cmath>-style free functions (ADL on MultiFloat)
// =============================================================================

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> abs(MultiFloat<T, N> const &x) {
  return std::signbit(x._limbs[0]) ? -x : x;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fabs(MultiFloat<T, N> const &x) {
  return abs(x);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fmin(MultiFloat<T, N> const &a,
                                MultiFloat<T, N> const &b) {
  return (a < b) ? a : b;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fmax(MultiFloat<T, N> const &a,
                                MultiFloat<T, N> const &b) {
  return (a < b) ? b : a;
}

template <typename T, std::size_t N>
constexpr bool signbit(MultiFloat<T, N> const &x) {
  return std::signbit(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr bool isfinite(MultiFloat<T, N> const &x) {
  for (std::size_t i = 0; i < N; ++i) {
    if (!std::isfinite(x._limbs[i])) {
      return false;
    }
  }
  return true;
}

template <typename T, std::size_t N>
constexpr bool isinf(MultiFloat<T, N> const &x) {
  return std::isinf(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr bool isnan(MultiFloat<T, N> const &x) {
  for (std::size_t i = 0; i < N; ++i) {
    if (std::isnan(x._limbs[i])) {
      return true;
    }
  }
  return false;
}

template <typename T, std::size_t N>
constexpr int fpclassify(MultiFloat<T, N> const &x) {
  return std::fpclassify(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> ldexp(MultiFloat<T, N> const &x, int n) {
  MultiFloat<T, N> r;
  for (std::size_t i = 0; i < N; ++i) {
    r._limbs[i] = std::ldexp(x._limbs[i], n);
  }
  return r;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> scalbn(MultiFloat<T, N> const &x, int n) {
  return ldexp(x, n);
}

template <typename T, std::size_t N>
constexpr int ilogb(MultiFloat<T, N> const &x) {
  return std::ilogb(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> copysign(MultiFloat<T, N> const &x,
                                    MultiFloat<T, N> const &y) {
  bool xs = std::signbit(x._limbs[0]);
  bool ys = std::signbit(y._limbs[0]);
  return (xs == ys) ? x : -x;
}

namespace detail {

// fast_two_sum (assumes |hi| >= |lo|), in-place renormalization helper.
template <typename T>
constexpr void renorm_fast(T &hi, T &lo) {
  T s = hi + lo;
  T b = s - hi;
  T e = lo - b;
  hi = s;
  lo = e;
}

} // namespace detail

// =============================================================================
// Rounding and integer-valued functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> floor(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::floor(x._limbs[0]);
  } else { // N == 2
    T fl_hi = std::floor(x._limbs[0]);
    if (fl_hi == x._limbs[0]) {
      // hi is already an integer; floor depends on the lo limb.
      r._limbs[0] = fl_hi;
      r._limbs[1] = std::floor(x._limbs[1]);
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      r._limbs[0] = fl_hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> ceil(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::ceil(x._limbs[0]);
  } else {
    T cl_hi = std::ceil(x._limbs[0]);
    if (cl_hi == x._limbs[0]) {
      r._limbs[0] = cl_hi;
      r._limbs[1] = std::ceil(x._limbs[1]);
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      r._limbs[0] = cl_hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> trunc(MultiFloat<T, N> const &x) {
  return std::signbit(x._limbs[0]) ? -floor(-x) : floor(x);
}

template <typename T, std::size_t N>
MultiFloat<T, N> round(MultiFloat<T, N> const &x) {
  // Round half away from zero, matching std::round.
  MultiFloat<T, N> half(T(0.5));
  return std::signbit(x._limbs[0]) ? -floor(-x + half) : floor(x + half);
}

template <typename T, std::size_t N>
MultiFloat<T, N> nearbyint(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::nearbyint(x._limbs[0]);
  } else {
    T hi = std::nearbyint(x._limbs[0]);
    if (hi == x._limbs[0]) {
      r._limbs[0] = hi;
      r._limbs[1] = std::nearbyint(x._limbs[1]);
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      r._limbs[0] = hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> rint(MultiFloat<T, N> const &x) {
  return nearbyint(x);
}

template <typename T, std::size_t N>
long lround(MultiFloat<T, N> const &x) {
  return std::lround(round(x)._limbs[0]);
}

template <typename T, std::size_t N>
long long llround(MultiFloat<T, N> const &x) {
  return std::llround(round(x)._limbs[0]);
}

template <typename T, std::size_t N>
long lrint(MultiFloat<T, N> const &x) {
  return std::lrint(rint(x)._limbs[0]);
}

template <typename T, std::size_t N>
long long llrint(MultiFloat<T, N> const &x) {
  return std::llrint(rint(x)._limbs[0]);
}

// =============================================================================
// Floating-point manipulation
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> frexp(MultiFloat<T, N> const &x, int *exp) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::frexp(x._limbs[0], exp);
  } else {
    int e;
    r._limbs[0] = std::frexp(x._limbs[0], &e);
    r._limbs[1] = std::ldexp(x._limbs[1], -e);
    *exp = e;
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> modf(MultiFloat<T, N> const &x, MultiFloat<T, N> *iptr) {
  *iptr = trunc(x);
  return x - *iptr;
}

template <typename T, std::size_t N>
MultiFloat<T, N> scalbln(MultiFloat<T, N> const &x, long n) {
  MultiFloat<T, N> r;
  for (std::size_t i = 0; i < N; ++i) {
    r._limbs[i] = std::scalbln(x._limbs[i], n);
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> logb(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  r._limbs[0] = std::logb(x._limbs[0]);
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> nextafter(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  if (x == y) {
    return y;
  }
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::nextafter(x._limbs[0], y._limbs[0]);
  } else {
    // Approximate one DD ulp ≈ ulp(hi) * 2^-53.
    T inf = std::numeric_limits<T>::infinity();
    T target = (x < y) ? inf : -inf;
    T next_hi = std::nextafter(x._limbs[0], target);
    T ulp = std::abs(next_hi - x._limbs[0]);
    T eps = std::ldexp(ulp, -53);
    return (x < y) ? x + MultiFloat<T, N>(eps) : x - MultiFloat<T, N>(eps);
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> nexttoward(MultiFloat<T, N> const &x,
                            MultiFloat<T, N> const &y) {
  return nextafter(x, y);
}

// =============================================================================
// Basic arithmetic helpers
// =============================================================================

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fma(MultiFloat<T, N> const &x,
                               MultiFloat<T, N> const &y,
                               MultiFloat<T, N> const &z) {
  // Not a hardware fma, but provides the cmath interface.
  return x * y + z;
}

template <typename T, std::size_t N>
MultiFloat<T, N> fmod(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  return x - trunc(x / y) * y;
}

template <typename T, std::size_t N>
MultiFloat<T, N> remainder(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  return x - round(x / y) * y;
}

template <typename T, std::size_t N>
MultiFloat<T, N> remquo(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y,
                        int *quo) {
  MultiFloat<T, N> q = round(x / y);
  *quo = static_cast<int>(q._limbs[0]);
  return x - q * y;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fdim(MultiFloat<T, N> const &x,
                                MultiFloat<T, N> const &y) {
  return (x > y) ? (x - y) : MultiFloat<T, N>();
}

// =============================================================================
// Power, exponential and logarithm
//
// For N == 2 these use a single first-order Newton/derivative correction
// against the leading-limb std:: implementation, which yields ~2x the
// precision of double but does not necessarily reach the full ~106-bit
// double-double precision.
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sqrt(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sqrt(x._limbs[0]);
    return r;
  } else {
    if (!(x._limbs[0] > T(0))) {
      r._limbs[0] = std::sqrt(x._limbs[0]); // 0, -0, NaN, or signaling
      return r;
    }
    // Karp/Markstein: r = s + (x - s*s) / (2s), evaluated in DD.
    T s = std::sqrt(x._limbs[0]);
    MultiFloat<T, N> s_dd(s);
    MultiFloat<T, N> residual = x - s_dd * s_dd;
    MultiFloat<T, N> correction(residual._limbs[0] * (T(0.5) / s));
    return s_dd + correction;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> cbrt(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::cbrt(x._limbs[0]);
    return r;
  } else {
    if (x._limbs[0] == T(0)) {
      return MultiFloat<T, N>();
    }
    T s = std::cbrt(x._limbs[0]);
    MultiFloat<T, N> s_dd(s);
    MultiFloat<T, N> residual = x - s_dd * s_dd * s_dd;
    MultiFloat<T, N> correction(residual._limbs[0] / (T(3) * s * s));
    return s_dd + correction;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> hypot(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  return sqrt(x * x + y * y);
}

template <typename T, std::size_t N>
MultiFloat<T, N> exp(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::exp(x._limbs[0]);
    return r;
  } else {
    // exp(hi + lo) = exp(hi) * (1 + lo + lo^2/2 + ...)
    T e = std::exp(x._limbs[0]);
    MultiFloat<T, N> e_dd(e);
    return e_dd + e_dd * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> exp2(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::exp2(x._limbs[0]);
    return r;
  } else {
    T e = std::exp2(x._limbs[0]);
    const T ln2 = std::log(T(2));
    MultiFloat<T, N> e_dd(e);
    // d/dx 2^x = 2^x * ln(2)
    return e_dd + e_dd * MultiFloat<T, N>(x._limbs[1] * ln2);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> expm1(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::expm1(x._limbs[0]);
    return r;
  } else {
    return exp(x) - MultiFloat<T, N>(T(1));
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log(x._limbs[0]);
    return r;
  } else {
    // log(hi + lo) = log(hi) + log(1 + lo/hi) ≈ log(hi) + lo/hi
    T l = std::log(x._limbs[0]);
    return MultiFloat<T, N>(l) +
           MultiFloat<T, N>(x._limbs[1] / x._limbs[0]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log10(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log10(x._limbs[0]);
    return r;
  } else {
    const T inv_ln10 = T(1) / std::log(T(10));
    return log(x) * MultiFloat<T, N>(inv_ln10);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log2(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log2(x._limbs[0]);
    return r;
  } else {
    const T inv_ln2 = T(1) / std::log(T(2));
    return log(x) * MultiFloat<T, N>(inv_ln2);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log1p(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log1p(x._limbs[0]);
    return r;
  } else {
    return log(x + MultiFloat<T, N>(T(1)));
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> pow(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::pow(x._limbs[0], y._limbs[0]);
    return r;
  } else {
    if (x._limbs[0] == T(0) && y._limbs[0] == T(0)) {
      return MultiFloat<T, N>(T(1));
    }
    return exp(y * log(x));
  }
}

// =============================================================================
// Trigonometric functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sin(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sin(x._limbs[0]);
    return r;
  } else {
    // sin(hi + lo) ≈ sin(hi) + cos(hi) * lo
    T s = std::sin(x._limbs[0]);
    T c = std::cos(x._limbs[0]);
    return MultiFloat<T, N>(s) +
           MultiFloat<T, N>(c) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> cos(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::cos(x._limbs[0]);
    return r;
  } else {
    T s = std::sin(x._limbs[0]);
    T c = std::cos(x._limbs[0]);
    return MultiFloat<T, N>(c) -
           MultiFloat<T, N>(s) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> tan(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::tan(x._limbs[0]);
    return r;
  } else {
    return sin(x) / cos(x);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> asin(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::asin(x._limbs[0]);
    return r;
  } else {
    // d/dx asin(x) = 1/sqrt(1 - x^2)
    T a = std::asin(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(MultiFloat<T, N>(T(1)) - x * x);
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> acos(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::acos(x._limbs[0]);
    return r;
  } else {
    // d/dx acos(x) = -1/sqrt(1 - x^2)
    T a = std::acos(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(MultiFloat<T, N>(T(1)) - x * x);
    return MultiFloat<T, N>(a) - MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> atan(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::atan(x._limbs[0]);
    return r;
  } else {
    // d/dx atan(x) = 1/(1 + x^2)
    T a = std::atan(x._limbs[0]);
    MultiFloat<T, N> denom = MultiFloat<T, N>(T(1)) + x * x;
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> atan2(MultiFloat<T, N> const &y, MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    return r;
  } else {
    // d atan2(y, x) = (x*dy - y*dx) / (x^2 + y^2)
    T a = std::atan2(y._limbs[0], x._limbs[0]);
    MultiFloat<T, N> num = x * MultiFloat<T, N>(y._limbs[1]) -
                           y * MultiFloat<T, N>(x._limbs[1]);
    MultiFloat<T, N> denom = x * x + y * y;
    return MultiFloat<T, N>(a) + num / denom;
  }
}

// =============================================================================
// Hyperbolic functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sinh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sinh(x._limbs[0]);
    return r;
  } else {
    T s = std::sinh(x._limbs[0]);
    T c = std::cosh(x._limbs[0]);
    return MultiFloat<T, N>(s) +
           MultiFloat<T, N>(c) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> cosh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::cosh(x._limbs[0]);
    return r;
  } else {
    T s = std::sinh(x._limbs[0]);
    T c = std::cosh(x._limbs[0]);
    return MultiFloat<T, N>(c) +
           MultiFloat<T, N>(s) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> tanh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::tanh(x._limbs[0]);
    return r;
  } else {
    return sinh(x) / cosh(x);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> asinh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::asinh(x._limbs[0]);
    return r;
  } else {
    // d/dx asinh(x) = 1/sqrt(1 + x^2)
    T a = std::asinh(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(MultiFloat<T, N>(T(1)) + x * x);
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> acosh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::acosh(x._limbs[0]);
    return r;
  } else {
    // d/dx acosh(x) = 1/sqrt(x^2 - 1)
    T a = std::acosh(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(x * x - MultiFloat<T, N>(T(1)));
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> atanh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::atanh(x._limbs[0]);
    return r;
  } else {
    // d/dx atanh(x) = 1/(1 - x^2)
    T a = std::atanh(x._limbs[0]);
    MultiFloat<T, N> denom = MultiFloat<T, N>(T(1)) - x * x;
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

// =============================================================================
// Error and gamma functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> erf(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::erf(x._limbs[0]);
    return r;
  } else {
    // d/dx erf(x) = 2/sqrt(pi) * exp(-x^2)
    const T two_over_sqrt_pi = T(2) / std::sqrt(std::acos(T(-1)));
    T e = std::erf(x._limbs[0]);
    T deriv = two_over_sqrt_pi * std::exp(-x._limbs[0] * x._limbs[0]);
    return MultiFloat<T, N>(e) +
           MultiFloat<T, N>(deriv) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> erfc(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::erfc(x._limbs[0]);
    return r;
  } else {
    const T two_over_sqrt_pi = T(2) / std::sqrt(std::acos(T(-1)));
    T e = std::erfc(x._limbs[0]);
    T deriv = -two_over_sqrt_pi * std::exp(-x._limbs[0] * x._limbs[0]);
    return MultiFloat<T, N>(e) +
           MultiFloat<T, N>(deriv) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> tgamma(MultiFloat<T, N> const &x) {
  // No std::digamma in <cmath>, so the N == 2 case falls back to a
  // leading-limb evaluation (single-double precision).
  MultiFloat<T, N> r;
  r._limbs[0] = std::tgamma(x._limbs[0]);
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> lgamma(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  r._limbs[0] = std::lgamma(x._limbs[0]);
  return r;
}

// =============================================================================
// Additional classification and ordered comparison
// =============================================================================

template <typename T, std::size_t N>
constexpr bool isnormal(MultiFloat<T, N> const &x) {
  return std::isnormal(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr bool isgreater(MultiFloat<T, N> const &x,
                         MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x > y);
}

template <typename T, std::size_t N>
constexpr bool isgreaterequal(MultiFloat<T, N> const &x,
                              MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x >= y);
}

template <typename T, std::size_t N>
constexpr bool isless(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x < y);
}

template <typename T, std::size_t N>
constexpr bool islessequal(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x <= y);
}

template <typename T, std::size_t N>
constexpr bool islessgreater(MultiFloat<T, N> const &x,
                             MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x != y);
}

template <typename T, std::size_t N>
constexpr bool isunordered(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  return isnan(x) || isnan(y);
}

using float64x2 = MultiFloat<double, 2>;

} // namespace multifloats
