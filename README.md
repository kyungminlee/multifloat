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

## Precision classes

The library distinguishes three classes of operations by the precision they
deliver:

| Class                                              | Typical max relative error | Operations |
| -------------------------------------------------- | -------------------------- | ---------- |
| **Full DD**                                        | ~1e-32 (one DD ulp)        | `+`, `-`, `*`, `/`, `sqrt`, `abs`, `neg`, `min`, `max`, `mod`, `modulo`, `dim`, `sign`, `hypot`, `aint`, `anint`, `fraction`, `scale`, `set_exponent`, `pow_int`, complex `+`/`-`/`*`, `cx_conjg`, `cx_abs`, `cx_aimag`, every constructor and assignment |
| **Single-double (derivative-corrected)**           | ~1e-16 (one double ulp)    | `exp`, `log`, `log10`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `erfc_scaled` |
| **Compound (chained derivative correction)**      | ~1e-12 to 1e-14            | `pow` (= `exp(b·log(a))`), `gamma`, `log_gamma`, `bessel_*`, complex transcendentals (`cx_exp`, `cx_log`, `cx_sin`, `cx_cos`, `cx_tan`, `cx_sinh`, `cx_cosh`, `cx_tanh`, `cx_asin`, `cx_acos`, `cx_atan`, `cx_asinh`, `cx_acosh`, `cx_atanh`), `cx_div`, `cx_sqrt` |

The transcendentals use a single first-order derivative correction over the
leading-limb `std::` (or libm) intrinsic. This roughly doubles the precision
of single-double for benign inputs but does not reach the full ~106-bit DD
precision of the arithmetic kernels.

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
