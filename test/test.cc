// Unit tests for multifloats.hh.
//
// Each public operation is tested over a curated list of inputs and the
// result is compared against a __float128 reference. A per-operation
// precision report is printed at the end.
//
// Built with g++ + libquadmath (Apple Clang lacks __float128).

#include "multifloats.hh"

#include <quadmath.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace mf = multifloats;
using MF2 = mf::MultiFloat<double, 2>;

// =============================================================================
// __float128 reference helpers
// =============================================================================

using q_t = __float128;

static q_t to_q(MF2 const &x) {
  return (q_t)x._limbs[0] + (q_t)x._limbs[1];
}

static MF2 from_q(q_t v) {
  double hi = (double)v;
  if (!std::isfinite(hi)) {
    MF2 r;
    r._limbs[0] = hi;
    r._limbs[1] = hi;
    return r;
  }
  double lo = (double)(v - (q_t)hi);
  // fast_two_sum normalization (|hi| >= |lo| holds by construction).
  double s = hi + lo;
  double err = lo - (s - hi);
  MF2 r;
  r._limbs[0] = s;
  r._limbs[1] = err;
  return r;
}

static double q_rel_err(q_t got, q_t expected) {
  if (expected == 0) {
    return (got == 0) ? 0.0 : (double)fabsq(got);
  }
  q_t diff = got - expected;
  if (diff < 0) {
    diff = -diff;
  }
  q_t mag = expected < 0 ? -expected : expected;
  return (double)(diff / mag);
}

static char const *qstr(q_t v) {
  static char buf[64];
  quadmath_snprintf(buf, sizeof(buf), "%.30Qg", v);
  return buf;
}

// =============================================================================
// Stats / reporting
// =============================================================================

struct Stats {
  double max_rel = 0.0;
  double sum_rel = 0.0;
  int count = 0;
  void update(double rel) {
    if (!std::isfinite(rel)) {
      return;
    }
    if (rel > max_rel) {
      max_rel = rel;
    }
    sum_rel += rel;
    ++count;
  }
};

static int g_failures = 0;
static int g_checks = 0;

#define REQUIRE(cond)                                                          \
  do {                                                                         \
    ++g_checks;                                                                \
    if (!(cond)) {                                                             \
      ++g_failures;                                                            \
      std::fprintf(stderr, "FAIL %s:%d  REQUIRE(%s)\n", __FILE__, __LINE__,    \
                   #cond);                                                     \
    }                                                                          \
  } while (0)

// Compare a multifloat result against a __float128 reference, accumulate
// statistics, and fail loudly if the relative error exceeds `tol`.
static void check_q(char const *op, MF2 const &got, q_t expected, double tol,
                    Stats &stats, MF2 const *a = nullptr,
                    MF2 const *b = nullptr) {
  ++g_checks;
  q_t got_q = to_q(got);
  double rel = q_rel_err(got_q, expected);
  stats.update(rel);
  if (!(rel <= tol)) {
    ++g_failures;
    std::fprintf(stderr, "FAIL [%s] rel_err=%g > tol=%g\n", op, rel, tol);
    if (a) {
      std::fprintf(stderr, "  a        = (%a, %a)\n", a->_limbs[0],
                   a->_limbs[1]);
    }
    if (b) {
      std::fprintf(stderr, "  b        = (%a, %a)\n", b->_limbs[0],
                   b->_limbs[1]);
    }
    std::fprintf(stderr, "  expected = %s\n", qstr(expected));
    std::fprintf(stderr, "  got      = %s\n", qstr(got_q));
    std::fprintf(stderr, "  got_limbs= (%a, %a)\n", got._limbs[0],
                 got._limbs[1]);
  }
}

// =============================================================================
// Tolerances
// =============================================================================
//
// MultiFloat<double, 2> has ~104 bits of precision. The kernels lose a few
// bits to ordering/normalization, so we accept up to ~2^-100 ≈ 8e-31 for
// linear ops and ~2^-95 for division (one Newton step from 1/y[0]).

static constexpr double kAddTol = 1e-30;
static constexpr double kMulTol = 1e-30;
static constexpr double kDivTol = 1e-28;

// =============================================================================
// Inputs
// =============================================================================

// All inputs are constructed by rounding a __float128 value down to a
// normalized double-double. This guarantees that to_q() is the exact
// inverse of the construction, so the __float128 reference is exact for
// every operation.
static std::vector<MF2> sample_inputs() {
  std::vector<MF2> v;
  auto push_q = [&v](q_t value) { v.push_back(from_q(value)); };

  // Pure doubles.
  push_q((q_t)0);
  push_q((q_t)1);
  push_q((q_t)-1);
  push_q((q_t)2);
  push_q((q_t)0.5);
  push_q((q_t)3);
  push_q((q_t)0.1);
  push_q((q_t)-0.1);
  push_q((q_t)123456789.0);

  // Irrational doubles with a non-trivial lo limb (1/3, sqrt(2), pi, e).
  push_q((q_t)1 / (q_t)3);
  push_q(sqrtq((q_t)2));
  push_q(M_PIq);
  push_q(-M_Eq);

  // 1 + tiny perturbations exercising the lo limb.
  push_q((q_t)1 + scalbnq((q_t)1, -53));
  push_q((q_t)1 - scalbnq((q_t)1, -54));
  push_q((q_t)1 + scalbnq((q_t)1, -80));

  // Large and small magnitudes (still representable as a normalized
  // double-double whose total span fits within float128 precision).
  push_q(scalbnq((q_t)1, 100));
  push_q(scalbnq((q_t)1, -100));
  push_q(scalbnq(M_PIq, 50));
  push_q(scalbnq(-M_Eq, -50));

  return v;
}

// =============================================================================
// Construction / conversion / comparison (exact properties)
// =============================================================================

static void test_construction_and_conversion() {
  MF2 z;
  REQUIRE(z._limbs[0] == 0.0);
  REQUIRE(z._limbs[1] == 0.0);

  MF2 one(1.0);
  REQUIRE(one._limbs[0] == 1.0);
  REQUIRE(one._limbs[1] == 0.0);
  REQUIRE(static_cast<double>(one) == 1.0);
  REQUIRE(to_q(one) == (q_t)1);

  MF2 negz(-0.0);
  REQUIRE(static_cast<double>(negz) == 0.0);
  REQUIRE(std::signbit(negz._limbs[0]));
}

static void test_equality_and_ordering() {
  for (MF2 const &a : sample_inputs()) {
    for (MF2 const &b : sample_inputs()) {
      q_t qa = to_q(a);
      q_t qb = to_q(b);
      REQUIRE((a == b) == (qa == qb));
      REQUIRE((a != b) == (qa != qb));
      REQUIRE((a < b) == (qa < qb));
      REQUIRE((a > b) == (qa > qb));
      REQUIRE((a <= b) == (qa <= qb));
      REQUIRE((a >= b) == (qa >= qb));
    }
  }
}

// =============================================================================
// Arithmetic — compared against __float128
// =============================================================================

static void test_unary(Stats &stats) {
  for (MF2 const &x : sample_inputs()) {
    q_t qx = to_q(x);
    check_q("neg", -x, -qx, 0.0, stats);
    check_q("pos", +x, qx, 0.0, stats);
  }
}

static void test_addition(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t expected = to_q(a) + to_q(b);
      check_q("add", a + b, expected, kAddTol, stats, &a, &b);

      MF2 c = a;
      c += b;
      check_q("add+=", c, expected, kAddTol, stats, &a, &b);
    }
  }
}

static void test_subtraction(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t expected = to_q(a) - to_q(b);
      check_q("sub", a - b, expected, kAddTol, stats);

      MF2 c = a;
      c -= b;
      check_q("sub-=", c, expected, kAddTol, stats);
    }
  }
}

static void test_multiplication(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t expected = to_q(a) * to_q(b);
      check_q("mul", a * b, expected, kMulTol, stats);

      MF2 c = a;
      c *= b;
      check_q("mul*=", c, expected, kMulTol, stats);
    }
  }
}

static void test_division(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      if (b._limbs[0] == 0.0) {
        continue;
      }
      q_t expected = to_q(a) / to_q(b);
      check_q("div", a / b, expected, kDivTol, stats);

      MF2 c = a;
      c /= b;
      check_q("div/=", c, expected, kDivTol, stats);
    }
  }
}

// =============================================================================
// cmath-style wrappers
// =============================================================================

static void test_abs_fmin_fmax(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &x : inputs) {
    q_t qx = to_q(x);
    q_t qexp = qx < 0 ? -qx : qx;
    check_q("abs", mf::abs(x), qexp, 0.0, stats);
    check_q("fabs", mf::fabs(x), qexp, 0.0, stats);
  }
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t qa = to_q(a);
      q_t qb = to_q(b);
      check_q("fmin", mf::fmin(a, b), qa < qb ? qa : qb, 0.0, stats);
      check_q("fmax", mf::fmax(a, b), qa < qb ? qb : qa, 0.0, stats);
    }
  }
}

static void test_copysign(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &x : inputs) {
    for (MF2 const &y : inputs) {
      q_t qx = to_q(x);
      bool xs = std::signbit(x._limbs[0]);
      bool ys = std::signbit(y._limbs[0]);
      q_t qexp = (xs == ys) ? qx : -qx;
      // Match the wrapper's exact bit-pattern semantics on the limb.
      check_q("copysign", mf::copysign(x, y), qexp, 0.0, stats);
    }
  }
}

static void test_classification() {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  MF2 finite(1.5);
  MF2 zero(0.0);
  MF2 neg(-1.5);
  MF2 pinf(inf), ninf(-inf), mnan(nan);

  REQUIRE(!mf::signbit(finite));
  REQUIRE(mf::signbit(neg));
  REQUIRE(mf::isfinite(finite));
  REQUIRE(!mf::isfinite(pinf));
  REQUIRE(!mf::isfinite(mnan));
  REQUIRE(mf::isinf(pinf));
  REQUIRE(mf::isinf(ninf));
  REQUIRE(!mf::isinf(finite));
  REQUIRE(mf::isnan(mnan));
  REQUIRE(!mf::isnan(finite));
  REQUIRE(mf::fpclassify(finite) == FP_NORMAL);
  REQUIRE(mf::fpclassify(zero) == FP_ZERO);
  REQUIRE(mf::fpclassify(pinf) == FP_INFINITE);
  REQUIRE(mf::fpclassify(mnan) == FP_NAN);

  // Non-canonical DDs where hi is (signed) zero: the value lives in lo,
  // so signbit/abs must consult lo rather than blindly trusting hi.
  MF2 tiny_neg;
  tiny_neg._limbs[0] = 0.0;
  tiny_neg._limbs[1] = -1e-300;
  REQUIRE(mf::signbit(tiny_neg));
  REQUIRE(!mf::signbit(mf::abs(tiny_neg)));

  MF2 tiny_pos;
  tiny_pos._limbs[0] = -0.0;
  tiny_pos._limbs[1] = 1e-300;
  REQUIRE(!mf::signbit(tiny_pos));
  REQUIRE(!mf::signbit(mf::abs(tiny_pos)));
}

static void test_ldexp_scalbn_ilogb(Stats &stats) {
  auto inputs = sample_inputs();
  int const ns[] = {-50, -1, 0, 1, 25};
  for (MF2 const &x : inputs) {
    if (!mf::isfinite(x) || x == MF2(0.0)) {
      continue;
    }
    q_t qx = to_q(x);
    for (int n : ns) {
      q_t qexp = qx * scalbnq((q_t)1, n);
      check_q("ldexp", mf::ldexp(x, n), qexp, 0.0, stats);
      check_q("scalbn", mf::scalbn(x, n), qexp, 0.0, stats);
    }
    // mf::ilogb is documented as delegating to the leading limb (matches
    // Julia/Fortran exponent semantics).
    REQUIRE(mf::ilogb(x) == std::ilogb(x._limbs[0]));
  }
}

// =============================================================================
// nextafter edge cases (power-of-2 boundaries)
// =============================================================================
//
// At |hi| = 2^k the downward ulp(hi) halves (2^(k-53) vs the upward 2^(k-52)),
// so the old formula eps = ldexp(ulp, -53) made the down-step 2× too small.
// The fix: always base eps on the UP-side ulp of |hi|, so the step is
// symmetric in both directions and nextafter(nextafter(x, +inf), -inf) == x.
static void test_nextafter_symmetry() {
  double const powers_of_two[] = {
      0x1p-50, 0x1p-10, 0x1p-1, 0x1p0, 0x1p1, 0x1p2, 0x1p3, 0x1p10, 0x1p50,
  };
  MF2 const pinf(std::numeric_limits<double>::infinity());
  MF2 const ninf(-std::numeric_limits<double>::infinity());

  for (double p2 : powers_of_two) {
    // Both signs — sign handling on the hi limb must not flip the ulp.
    for (double sign : {+1.0, -1.0}) {
      MF2 x(sign * p2);
      MF2 up = mf::nextafter(x, pinf);
      MF2 dn = mf::nextafter(x, ninf);

      // Both steps must be nonzero (no collapse to x).
      REQUIRE(up != x);
      REQUIRE(dn != x);
      // Step directions are correct.
      REQUIRE(up > x);
      REQUIRE(dn < x);
      // Round-trip identity: step and return lands exactly on x.
      REQUIRE(mf::nextafter(up, ninf) == x);
      REQUIRE(mf::nextafter(dn, pinf) == x);
      // Step magnitudes are symmetric in DD (not halved on one side).
      q_t qx = to_q(x);
      q_t qup = to_q(up);
      q_t qdn = to_q(dn);
      q_t up_step = qup - qx;
      q_t dn_step = qx - qdn;
      // Equality as q_t: we set up_step/dn_step to match bit-for-bit.
      REQUIRE(up_step == dn_step);
      // Step equals expected eps = 2^-105 * p2 (i.e., DD-rel ulp × |x|).
      q_t expected = scalbnq((q_t)p2, -105);
      REQUIRE(up_step == expected);
    }
  }

  // Non-power-of-2 values: step up and down should be equal (trivially —
  // ulp is unambiguous there) and round-trip still identity.
  double const non_p2[] = {
      1.5, 1.25, 1.125, 3.0, 5.0, 7.0, 0.75, 0.375, 0x1.23456789abcdep10,
  };
  for (double v : non_p2) {
    for (double sign : {+1.0, -1.0}) {
      MF2 x(sign * v);
      MF2 up = mf::nextafter(x, pinf);
      MF2 dn = mf::nextafter(x, ninf);
      REQUIRE(up != x);
      REQUIRE(dn != x);
      REQUIRE(mf::nextafter(up, ninf) == x);
      REQUIRE(mf::nextafter(dn, pinf) == x);
      q_t up_step = to_q(up) - to_q(x);
      q_t dn_step = to_q(x) - to_q(dn);
      REQUIRE(up_step == dn_step);
    }
  }

  // Identity: x == y implies nextafter returns y unchanged (exact).
  for (double p2 : powers_of_two) {
    MF2 x(p2);
    REQUIRE(mf::nextafter(x, x) == x);
  }

  // Non-trivial DD inputs (hi and lo both present): round-trip should hold.
  q_t const seeds[] = {
      (q_t)1 + scalbnq((q_t)1, -60),
      M_PIq,
      sqrtq((q_t)2) * scalbnq((q_t)1, 30),
      (q_t)1 / (q_t)7,
  };
  for (q_t seed : seeds) {
    MF2 x = from_q(seed);
    MF2 up = mf::nextafter(x, pinf);
    MF2 dn = mf::nextafter(x, ninf);
    REQUIRE(mf::nextafter(up, ninf) == x);
    REQUIRE(mf::nextafter(dn, pinf) == x);
  }
}

// =============================================================================
// atan cutover edge cases
// =============================================================================
//
// The large-argument branch uses `t = -1/|x|` once `|x| >= 32/3`, which puts
// worst-case `|t| = 3/32 = 0.09375` at the fitted edge of the rational fit.
// The table branch spans k = 0..85 so the newly-reached cells (k = 83, 84, 85)
// and all crossings between them need to be walked.
static void test_atan_cutover(Stats &stats) {
  // Enumerate exact boundary values (k/8 centers where t = 0 exactly),
  // midpoints between them (worst table-branch |t|), and the new cutover.
  double const boundaries[] = {
      // Old cutover, now in table branch (k = 82 center).
      10.25,
      // New table entries: k/8 centers.
      10.375,                 // k = 83
      10.5,                   // k = 84
      10.625,                 // k = 85
      // Table-branch midpoints (maximum |t| within a cell).
      10.3125,                // between k=82 and k=83
      10.4375,                // between k=83 and k=84
      10.5625,                // between k=84 and k=85
      // New cutover: |x| = 32/3 -> large branch, |t| = 3/32 = 0.09375 exactly.
      32.0 / 3.0,
      // Large-branch values immediately beyond the cutover.
      10.7, 11.0, 12.0, 15.0, 100.0, 1e10, 1e100, 1e300,
      // Also pick up a value in the old 10.25..10.6667 gap that used to go
      // through the large branch at |t| > 0.09375.
      10.26, 10.3, 10.4, 10.5, 10.6, 10.65,
  };

  auto check_atan = [&](q_t qx) {
    MF2 x = from_q(qx);
    MF2 got = mf::atan(x);
    q_t expected = atanq(qx);
    // atan is always well-conditioned (|d/dx atan| <= 1), so tolerance is
    // the standard DD kernel tier.
    check_q("atan_cutover", got, expected, 1e-30, stats, &x);
  };

  for (double b : boundaries) {
    // The point itself, its bit-neighbors, and a ±16 ulp ring on the DD-side.
    q_t qb = (q_t)b;
    check_atan(qb);
    check_atan(-qb);
    q_t up = qb, dn = qb;
    for (int i = 0; i < 16; ++i) {
      up = nextafterq(up, (q_t)std::numeric_limits<double>::infinity());
      dn = nextafterq(dn, (q_t)-std::numeric_limits<double>::infinity());
      check_atan(up);
      check_atan(dn);
      check_atan(-up);
      check_atan(-dn);
    }
  }

  // Dense sweep across the newly-reached table range, with nontrivial lo
  // limbs so the kernel sees full-DD arguments (not just exact doubles).
  for (int i = 0; i < 2000; ++i) {
    // Uniform random in [10.25, 32/3).
    double frac = (double)i / 2000.0;
    q_t qx = (q_t)10.25 + ((q_t)(32.0 / 3.0) - (q_t)10.25) * (q_t)frac;
    // Add a sub-ulp perturbation so the lo limb is exercised.
    qx += scalbnq((q_t)1, -60) * (q_t)((i % 7) - 3);
    check_atan(qx);
    check_atan(-qx);
  }

  // Verify that the cutover value lands in the large branch: round-trip
  // identity arctan(x) + arctan(1/x) = pi/2 for x > 0 at |x| = 32/3.
  MF2 x = from_q((q_t)(32.0 / 3.0));
  MF2 one(1.0);
  MF2 y = one / x;
  MF2 sum = mf::atan(x) + mf::atan(y);
  q_t got = to_q(sum);
  q_t expected = M_PIq / (q_t)2;
  q_t diff = got - expected;
  if (diff < 0) diff = -diff;
  q_t mag = expected < 0 ? -expected : expected;
  double rel = (double)(diff / mag);
  ++g_checks;
  stats.update(rel);
  if (!(rel <= 1e-30)) {
    ++g_failures;
    std::fprintf(stderr,
                 "FAIL [atan+arccotan=pi/2 at 32/3] rel_err=%g\n", rel);
  }
}

// =============================================================================
// Driver
// =============================================================================

static void print_stats(char const *name, Stats const &s) {
  if (s.count == 0) {
    std::printf("  %-12s n=%d\n", name, s.count);
    return;
  }
  std::printf("  %-12s n=%-6d max_rel=%.3e  mean_rel=%.3e\n", name, s.count,
              s.max_rel, s.sum_rel / s.count);
}

int main() {
  test_construction_and_conversion();
  test_equality_and_ordering();
  test_classification();
  test_nextafter_symmetry();

  Stats add, sub, mul, div, una, abs_fmm, csgn, ldx, atn;
  test_unary(una);
  test_addition(add);
  test_subtraction(sub);
  test_multiplication(mul);
  test_division(div);
  test_abs_fmin_fmax(abs_fmm);
  test_copysign(csgn);
  test_ldexp_scalbn_ilogb(ldx);
  test_atan_cutover(atn);

  std::printf("[multifloats_test] %d checks, %d failures\n", g_checks,
              g_failures);
  std::printf("[multifloats_test] precision report (relative error vs "
              "__float128):\n");
  print_stats("unary", una);
  print_stats("add", add);
  print_stats("sub", sub);
  print_stats("mul", mul);
  print_stats("div", div);
  print_stats("abs/fmin/fmax", abs_fmm);
  print_stats("copysign", csgn);
  print_stats("ldexp/etc", ldx);
  print_stats("atan cutover", atn);

  return g_failures == 0 ? 0 : 1;
}
