// Decimal I/O for float64x2 — scientific-notation formatting.
//
// Layering. The core formatter lives as `multifloats::detail::format_scientific_chars`
// — a named C++ function, called directly by the header-inline
// `multifloats::to_chars` in multifloats.hh. The extern "C" `to_charsdd`
// entry point (declared in multifloats_c.h) is a one-line shim on top of
// the same detail primitive, for C / Fortran consumers. This avoids the
// C++ → C → C++ round-trip that a straight extern-"C"-as-the-primary entry
// would impose on header-inline C++ callers.
//
// The `multifloats::to_string` wrapper (also inline in multifloats.hh) sits
// above `to_chars`, so the returned `std::string` is constructed entirely
// in the consumer's translation unit — keeping the library's link surface
// free of `std::string` and the libstdc++ `_GLIBCXX_USE_CXX11_ABI` dual-ABI
// mangling distinction.

#include "multifloats.hh"

#include <cmath>
#include <cstddef>
#include <cstdio>   // std::snprintf for exponent formatting
#include <cstring>
#include <ostream>

namespace {

// Renormalize (hi, lo) into canonical DD form.
inline void io_renorm(double &hi, double &lo) {
  double s = hi + lo;
  double err = lo - (s - hi);
  hi = s;
  lo = err;
}

// (hi, lo) *= d, with d a double; result renormalized.
inline void io_dd_mul_d(double &hi, double &lo, double d) {
  double p = hi * d;
  double e = std::fma(hi, d, -p);
  lo = lo * d + e;
  hi = p;
  io_renorm(hi, lo);
}

} // anonymous namespace

namespace multifloats {
namespace detail {

// Scientific-notation decimal formatter for float64x2. Writes into
// [first, last) without a NUL terminator, returns one-past-last-written
// on success or nullptr when the output wouldn't fit (in which case the
// buffer is left unchanged). `precision` is clamped to [1, 34]. Special
// values emit "nan", "inf", "-inf", "0e+00", "-0e+00".
//
// Strategy: produce the full output into a fixed 48-byte tmp[] buffer
// (sized to hold any legal output), then size-check once and memcpy into
// the caller's range. Any legal output fits in at most ~42 bytes.
char *format_scientific_chars(float64x2 const &value, int precision,
                              char *first, char *last) noexcept {
  double hi = value._limbs[0];
  double lo = value._limbs[1];

  char tmp[48];
  std::size_t w = 0;

  auto put_char = [&](char c) {
    if (w < sizeof(tmp)) tmp[w++] = c;
  };
  auto put_str = [&](const char *s) {
    while (*s && w < sizeof(tmp)) tmp[w++] = *s++;
  };

  if (std::isnan(hi) || std::isnan(lo)) {
    put_str("nan");
  } else if (std::isinf(hi)) {
    put_str(hi > 0 ? "inf" : "-inf");
  } else if (hi == 0.0) {
    put_str(std::signbit(hi) ? "-0e+00" : "0e+00");
  } else {
    if (precision < 1) precision = 1;
    if (precision > 34) precision = 34;

    const bool neg = hi < 0.0;
    if (neg) { hi = -hi; lo = -lo; }

    // Decimal exponent estimate. log10 may be off by 1 due to rounding; a
    // post-scale fixup corrects.
    int e10 = static_cast<int>(std::floor(std::log10(hi)));

    // Scale to [1, 10). Multiplying by 10 is exact in DD; multiplying by 0.1
    // is not, but each iteration retains ~106 bits of precision, so 308
    // iterations (max for dp range) lose at most a handful of ulps — well
    // below the visible 34 digits.
    int shift = -e10;
    while (shift > 0) { io_dd_mul_d(hi, lo, 10.0); --shift; }
    while (shift < 0) { io_dd_mul_d(hi, lo, 0.1);  ++shift; }
    // Drift correction: at most ±1 from the log10 estimate.
    if (hi >= 10.0)     { io_dd_mul_d(hi, lo, 0.1);  ++e10; }
    else if (hi < 1.0)  { io_dd_mul_d(hi, lo, 10.0); --e10; }

    // Extract precision+2 digits; the last two guard round-half-to-even.
    const int ndigits = precision + 2;
    char digits[36] = {};
    for (int i = 0; i < ndigits; ++i) {
      int d = static_cast<int>(hi);
      if (d < 0) d = 0;
      else if (d > 9) d = 9;
      digits[i] = static_cast<char>('0' + d);
      hi -= d;
      io_renorm(hi, lo);
      io_dd_mul_d(hi, lo, 10.0);
    }

    int guard = (digits[precision] - '0') * 10 + (digits[precision + 1] - '0');
    bool round_up = guard > 50 ||
                    (guard == 50 && ((digits[precision - 1] - '0') & 1));
    if (round_up) {
      for (int i = precision - 1; i >= 0; --i) {
        if (digits[i] < '9') { ++digits[i]; break; }
        digits[i] = '0';
        if (i == 0) {
          std::memmove(digits + 1, digits,
                       static_cast<std::size_t>(precision) - 1);
          digits[0] = '1';
          ++e10;
        }
      }
    }
    digits[precision] = '\0';

    if (neg) put_char('-');
    put_char(digits[0]);
    if (precision > 1) {
      put_char('.');
      for (int i = 1; i < precision; ++i) put_char(digits[i]);
    }
    put_char('e');
    int abs_e = e10 < 0 ? -e10 : e10;
    put_char(e10 < 0 ? '-' : '+');
    // Two-digit minimum exponent (matches libc `%e` formatting).
    char ebuf[8];
    int ew = std::snprintf(ebuf, sizeof(ebuf), "%d", abs_e);
    if (ew == 1) put_char('0');
    for (int i = 0; i < ew; ++i) put_char(ebuf[i]);
  }

  std::size_t avail = (last > first) ? static_cast<std::size_t>(last - first) : 0;
  if (w > avail) return nullptr;
  std::memcpy(first, tmp, w);
  return first + w;
}

} // namespace detail
} // namespace multifloats

extern "C" {
char *to_charsdd(float64x2_t x, int precision, char *first, char *last) {
  return multifloats::detail::format_scientific_chars(
      multifloats::detail::from_f64x2(x), precision, first, last);
}
} // extern "C"

namespace multifloats {

std::ostream &operator<<(std::ostream &os, float64x2 const &x) {
  int p = static_cast<int>(os.precision());
  if (p <= 17) p = 32;
  char buf[MULTIFLOATS_DD_CHARS_BUFSIZE];
  char *end = detail::format_scientific_chars(x, p, buf, buf + sizeof(buf));
  // `buf` is sized to hold any legal output, so `end` is non-null here.
  if (end) os.write(buf, end - buf);
  return os;
}

} // namespace multifloats
