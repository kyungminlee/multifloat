# multifloats

Double-double arithmetic for Fortran and C++.

`multifloats` provides two derived/template types,

| Type             | Storage             | Precision        |
| ---------------- | ------------------- | ---------------- |
| `float64x2`      | 2 × `real(dp)`      | ~106 bits (~32 decimal digits) |
| `complex128x2`   | 2 × `float64x2`     | ~106 bits per component |

implemented natively from error-free transformations on regular IEEE
doubles. There are no quad-precision (`real(16)` / `__float128`) temporaries
in any arithmetic, transcendental, or array kernel — those are reserved
exclusively for the `to_qp` / `from_qp` conversion helpers and the defined
I/O routines.

The Fortran module ships in `fsrc/multifloats.fypp` (preprocessed via
[fypp]) and the C++ header in `src/multifloats.hh`. Both expose the same
algorithmic surface so the same kernels can be used from either language.

## Features

### Fortran (`fsrc/multifloats.fypp` → `multifloats` module)

- **Types** — `float64x2` and `complex128x2`. Both carry the `SEQUENCE`
  attribute so they can appear in `EQUIVALENCE` statements (required by
  some LAPACK routines such as `wlaln2` / `DLALN2`).
- **Operators** — `+`, `-`, `*`, `/`, `**`, `==`, `/=`, `<`, `>`, `<=`,
  `>=` between every combination of {`float64x2`, `complex128x2`,
  `real(dp)`, `real(sp)`, `integer`, `complex(dp)`, `complex(sp)`}.
- **Constructors and assignment** — every supported numeric kind, in
  both directions, including the non-default integer kinds
  `integer(int8)`, `integer(int16)`, `integer(int64)`. Identity
  constructors `float64x2(float64x2)` and `complex128x2(complex128x2)`
  let generic code call the constructor regardless of input type.
- **`<cmath>`-style intrinsics** — `abs`, `sqrt`, `cbrt`, `exp`, `exp2`,
  `expm1`, `log`, `log2`, `log10`, `log1p`, `sin`, `cos`, `tan`, `asin`,
  `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`,
  `atanh`, `erf`, `erfc`, `erfc_scaled`, `gamma`, `log_gamma`, `pow`,
  `hypot`, `fmod`/`mod`, `modulo`, `dim`, `sign`, `min`, `max` (3..8
  arguments), `bessel_j0/j1/y0/y1/jn/yn`, `aint`, `anint`, `floor`,
  `ceiling`, `nint`, `int`, `dble`, `fraction`, `spacing`, `rrspacing`,
  `scale`, `set_exponent`, `exponent`, `nearest`, `epsilon`, `huge`,
  `tiny`, `precision`, `range`, `digits`, `radix`, `minexponent`,
  `maxexponent`, `storage_size`.
- **Reductions** — `sum`, `product`, `maxval`, `minval`, `maxloc`,
  `minloc`, `findloc`, `dot_product`, `norm2`, `matmul` for both
  `float64x2` and `complex128x2`, supporting all rank-1..7 forms with
  `dim=`, `mask=`, and `back=` arguments.
- **Other** — `random_number` (rank 0..7), defined formatted I/O for
  both types, and inquiry functions matching `real(dp)` semantics on the
  leading limb.

The `float64x2` interface is designed to mirror `REAL(KIND=16)` so that
existing quad-precision code can switch types with minimal changes.

### C++ (`src/multifloats.hh`)

A single C++17 header in the `multifloats` namespace, providing
`MultiFloat<T, N>` (`N` ∈ {1, 2}) with the same algorithmic kernels and a
convenience alias `using float64x2 = MultiFloat<double, 2>;`.

The full `<cmath>` double-name surface (`floor`, `ceil`, `trunc`, `round`,
`nearbyint`, `rint`, `lround`, `llround`, `lrint`, `llrint`, `frexp`,
`modf`, `scalbn`, `scalbln`, `ilogb`, `logb`, `nextafter`, `nexttoward`,
`copysign`, `fma`, `fmod`, `remainder`, `remquo`, `fdim`, `sqrt`, `cbrt`,
`hypot`, `pow`, `exp`, `exp2`, `expm1`, `log`, `log10`, `log2`, `log1p`,
`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`,
`tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `tgamma`, `lgamma`,
`abs`, `fabs`, `fmin`, `fmax`, `fpclassify`, `isfinite`, `isinf`, `isnan`,
`isnormal`, `signbit`, `isgreater`, `isgreaterequal`, `isless`,
`islessequal`, `islessgreater`, `isunordered`) is implemented via ADL on
`MultiFloat`. The C++ side has zero external dependencies — `__float128`
is only used by the test harness as a high-precision reference.

## Precision

`multifloat_fuzz` runs every operation listed below through 1M random
input pairs (using adversarial input strategies that include subnormals,
near-cancellation pairs, near-overflow pairs, and non-finite leading limbs)
and reports the per-op `(max_rel_err, mean_rel_err)` against a
quad-precision (`real(16)`) reference. The numbers in the tables below are
representative samples from a recent run of `multifloat_fuzz` (seed = 0,
1M iterations); your build will land in the same orders of magnitude.

A relative error reported as `0` means *exactly bit-equal* to the qp
reference for every input the fuzz drew. Operations are organized by their
precision class.

### Full double-double (~1e-32 max — one DD ulp)

The arithmetic kernels and operations whose result is a bit-exact rearrangement
of the input limbs. Worst-case error is one DD ulp, mean error around 1/100
of a DD ulp.

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `+` (mf±mf, mf±dp, etc.) | 1.5e-32 | 3.0e-34 |
| `-` | 6.2e-33 | 1.8e-34 |
| `*` | 3.5e-32 | 6.5e-34 |
| `/` | 5.3e-32 | 1.8e-33 |
| `sqrt` (Karp/Markstein iteration) | 5.0e-32 | 5.3e-33 |
| `mod`, `modulo` | 2.2e-32 | 9e-36 |
| `dim` | 6.2e-33 | 9.4e-35 |
| `min`, `max`, `min(a..h)`, `max(a..h)` | 6.2e-33 | 9.5e-35 |
| `hypot` (with overflow-safe scaling) | 7.8e-32 | 5.2e-33 |
| `pow` (integer exponent, repeated multiplication) | 2.2e-32 | 1.3e-33 |
| Complex `+`, `-` (real and imag parts) | 1.4e-32 | 3.2e-34 |
| Complex `*` real part | 0 | 0 |
| Complex `*` imag part | 1.9e-32 | 1.5e-33 |
| Complex `/` real part | 4.5e-32 | 2.5e-33 |
| `cx_conjg`, `cx_abs`, `cx_aimag` | 6.7e-32 | 5.4e-33 |

### Bit-exact (always 0 error)

Operations that are pure limb manipulation, sign flips, or trivial
promotions / truncations. The fuzz reports `max_rel = mean_rel = 0` over
1M iterations for every operation in this group.

- **Unary**: `abs`, `neg`, `sign`, `aint`, `anint`, `fraction`,
  `scale`, `set_exponent`
- **Mixed-mode arithmetic where one side is dp**: `mf + dp` (`add_fd`)
  reports 0 because the lo-limb error is in the dp ulp range
- **Every constructor**: `float64x2(...)` and `complex128x2(...)` for
  every supported numeric kind
- **Every assignment**: `mf ↔ {dp, sp, int, int8, int16, int64, cdp, csp}`
  and `cx ↔ {dp, sp, int, int8, int16, int64, cdp, csp}`
- **Complex `*` real part**: bit-exact because the real part is computed
  as a single fma-style chain with no cancellation between terms

### Single-double, first-order derivative corrected (~1e-16 max)

Real transcendentals computed as `f(hi) + f'(hi) · lo` combined via
`fast_two_sum`. Each gives roughly one double ulp of relative error,
which is single-double precision — *not* full DD precision. The cost
saving (no quad-precision temporaries, no polynomial tables) is the
deliberate tradeoff.

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `exp` | 1.1e-16 | 2.2e-17 |
| `log` | 1.1e-16 | 2.4e-17 |
| `log10` | 1.0e-16 | 2.4e-17 |
| `sin` | 1.1e-16 | 2.3e-17 |
| `cos` | 1.1e-16 | 1.8e-17 |
| `tan` | 2.9e-16 | 3.8e-17 |
| `asin` | 1.2e-16 | 1.8e-17 |
| `acos` | 1.6e-16 | 4.4e-17 |
| `atan` | 1.3e-16 | 1.6e-17 |
| `atan2` | 1.0e-16 | 1.4e-17 |
| `sinh` | 1.1e-16 | 2.1e-17 |
| `cosh` | 1.1e-16 | 2.7e-17 |
| `tanh` | 1.1e-16 | 1.8e-17 |
| `asinh` | 1.9e-16 | 2.2e-17 |
| `acosh` | 1.7e-16 | 1.6e-17 |
| `atanh` | 2.1e-16 | 2.5e-17 |
| `erf` | 3.8e-16 | 1.0e-16 |
| `erfc` | 5.2e-16 | 3.6e-17 |
| `erfc_scaled` | 4.1e-16 | 4.8e-17 |

### Compound — chained derivative corrections (~1e-12 to 1e-14 max)

Functions that internally chain two or more single-double-precision
transcendentals or have cancellation in their construction. Worst-case
error is several orders of magnitude looser than a single transcendental,
but mean precision is typically still ~1e-15.

| Op | max_rel | mean_rel | Notes |
| --- | --- | --- | --- |
| `pow` (mf**mf), `mf**dp`, `dp**mf` | 3.4e-14 | 5.7e-16 | implemented as `exp(b·log(a))` |
| `gamma` | 1.2e-16 | 4.0e-17 | derivative correction unavailable; uses libm directly |
| `log_gamma` | 3.3e-16 | 5.0e-17 | likewise |
| `bessel_j0` | 3.9e-16 | 2.4e-17 | libm Bessel precision |
| `bessel_j1` | 2.5e-15 | 3.1e-17 | |
| `bessel_jn` (n=3 sample) | 5.0e-15 | 4.4e-17 | |
| `bessel_y0` | 2.3e-15 | 8.3e-17 | |
| `bessel_y1` | 7.9e-16 | 8.3e-17 | |
| `bessel_yn` (n=3 sample) | 1.4e-13 | 2.1e-16 | |
| `cx_exp` | 1.9e-16 | 5.0e-17 | derived from real `exp`, `sin`, `cos` |
| `cx_log` | 1.1e-16 | 2.0e-17 | |
| `cx_sin`, `cx_cos`, `cx_sinh`, `cx_cosh` | 2.2e-16 | 5.1e-17 | |
| `cx_tan`, `cx_tanh` | 1.9e-14 | 1.4e-16 | catastrophic cancellation in `sin/cos` near poles |
| `cx_atan` | 1.1e-16 | 3.9e-17 | |
| `cx_asin`, `cx_acos`, `cx_atanh`, `cx_acosh` | 1.1e-16 | 3.9e-17 | range-restricted in fuzz to avoid the near-zero cancellation regime |
| `cx_asinh` | 1.2e-16 | 3.9e-17 | |
| `cx_div` real part | 4.5e-32 | 2.5e-33 | full DD on real part |
| `cx_div` imag part | 8.7e-16 | 6.8e-19 | cancellation |
| `cx_sqrt` real / imag (with even-power scaling) | 2.5e-20 | 1e-23 | almost full DD; loses a few digits to the cancellation in `sqrt((|z|±a)/2)` |

### Array reductions (small-array fuzz, n = 8)

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `sum`     | 4.0e-30 | 1.3e-32 |
| `product` | 1.1e-49 | 2.6e-52 |
| `maxval`  | 5.9e-33 | 1.5e-33 |
| `minval`  | 6.1e-33 | 1.4e-33 |
| `dot_product` | 1.2e-31 | 7.0e-33 |
| `norm2`   | 6.0e-32 | 1.0e-32 |
| `matmul` (mv, n=8) | 5.0e-31 | 6.6e-33 |

The reductions accumulate over `n` elements, so worst-case error grows
linearly with `n` while staying inside the full-DD regime.

### Why some functions only get single-double precision

The transcendentals (`exp`, `log`, `sin`, `cos`, ...) are evaluated as

  `f(hi + lo) ≈ f(hi) + f'(hi) · lo`

where `f(hi)` and `f'(hi)` come from the standard library on the leading
limb only. This is *one* Newton step rather than the full polynomial
evaluation that
[Julia's MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl)
uses to reach full DD precision in `exp`/`log`. The trade-off:

- **No quad-precision dependency** — there are no polynomial coefficient
  tables and no `__float128` / libquadmath calls in any kernel.
- **Roughly doubles the precision of double** for benign inputs, but
  cancellation in the formula or in `f'(hi)·lo` itself can reduce the
  effective precision back toward single-double in pathological cases.

If you need true ~106-bit transcendentals, the route is to port the Julia
polynomial evaluators (the tables under `external/MultiFloats.jl/src/`)
into the `mf_${func}$` cases of `fsrc/multifloats.fypp`. The rest of the
infrastructure (operators, reductions, complex / array support) is
already at full DD.

## Building

Requires:
- CMake ≥ 3.27
- A Fortran 2018 compiler with `REAL(KIND=16)` support and a C++17
  compiler with `__float128` / libquadmath. On macOS this means Homebrew
  GCC 13/14/15 (Apple Clang and Apple-shipped LLVM Flang are not
  sufficient). The build pins `g++-15` / `gfortran-15` automatically; pass
  `-DCMAKE_CXX_COMPILER=...` / `-DCMAKE_Fortran_COMPILER=...` to override.
- [`fypp`](https://fypp.readthedocs.io/) on `PATH` (the Fortran source is
  generated at configure time).

```sh
cmake -B build -S .
cmake --build build
ctest --test-dir build --output-on-failure
```

The build produces:
- `libmultifloat-fortran.a` — the Fortran module library (and the
  generated `.mod` files in `build/modules/`).
- `libblas-multifloat.a` — `wgemm` / `wtrsm` BLAS shims that operate on
  `float64x2` matrices.
- Five test executables (see below).

## Tests

```sh
ctest --test-dir build --output-on-failure
```

| Test                  | Language | What it covers |
| --------------------- | -------- | -------------- |
| `multifloats_precision` | Fortran  | Targeted vs-quad precision checks for constructors, assignments, every arithmetic and reduction op, and edge-case sweeps (signed zero, infinities, NaN propagation, subnormal/huge boundary, ULP boundary, dim/mask/back reduction variants). |
| `multifloat_fuzz`       | Fortran  | 1M random pairs through every public function for which random real input is meaningful. Adversarial input strategies cover subnormals, near-cancellation, overflow boundary, and non-finite limbs. Prints a per-op `(max_rel, mean_rel, count)` precision report at the end. |
| `multifloat_test`       | Fortran  | Hand-written assertions for arithmetic, signed zero, NaN/Inf propagation, classification, math intrinsics, and rounding. |
| `multifloats_test`      | C++      | Targeted vs-`__float128` precision checks for `src/multifloats.hh`. |
| `multifloats_fuzz`      | C++      | C++ fuzz with the same precision-report machinery as the Fortran fuzz. |

## Layout

```
fsrc/multifloats.fypp     -- Fortran source (fypp template)
src/multifloats.hh        -- C++17 header
blas/                     -- BLAS shims for float64x2
test/                     -- Fortran and C++ test suites
external/                 -- vendored references (MultiFloats.jl, LAPACK)
```

The C++ kernels follow the algorithms in
[Julia's MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl) (vendored
under `external/`); the Fortran kernels are direct translations of the same
double-double algorithms.

## License

See [LICENSE](LICENSE) (if present) or the source headers for licensing terms.

[fypp]: https://fypp.readthedocs.io/
