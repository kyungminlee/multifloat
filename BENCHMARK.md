# Benchmark Results

Comparison of multifloats double-double (DD) arithmetic against quad
precision (`real(16)` / `__float128` via libquadmath). The multifloats
implementations use explicit two-limb error-free transformations (EFTs)
on pairs of `double` values, achieving ~106 bits of significand — close
to quad precision's 113 bits — at a fraction of the cost.

## System

| | |
|---|---|
| **CPU** | Apple M1 Max (ARM64, 10 cores) |
| **RAM** | 64 GB |
| **OS** | macOS 26.3 (Darwin 25.3.0) |
| **Fortran** | GNU Fortran (Homebrew GCC 15.2.0\_1) 15.2.0 |
| **C++** | g++-15 (Homebrew GCC 15.2.0\_1) 15.2.0 |
| **Build** | CMake 4.3.1, `-O3 -flto`, OBJECT library (not STATIC — macOS `ar` strips GIMPLE from Mach-O, breaking LTO through `.a` archives) |

## Precision key

Precision is measured as the maximum relative error vs the quad-precision
(`real(16)`) reference over ~1M random inputs. The "precision" column
indicates the kernel implementation strategy:

| Label | max\_rel | Meaning |
|---|---|---|
| **full DD** | ~1e-30 | Full double-double EFT kernel (~106 bits) |
| **exact** | 0.0 | Bit-exact (no rounding involved) |
| **deriv-corrected** | ~1e-25 to 1e-18 | `f(hi) + f'(hi)*lo` correction gives near-DD |
| **single-double** | ~1e-16 to 1e-14 | Leading-limb libm call, no lo correction |

Operations marked **single-double** trade DD precision for massive speed
gains (60–165×). They use the leading `dp` limb only for the function
evaluation, achieving `double`-level accuracy (~15 digits) rather than
full DD (~31 digits).

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. The "speedup" column is `qp_time / mf_time`:
values > 1× mean multifloats is faster.

### Arithmetic

| op | speedup | max\_rel | precision |
|---|---|---|---|
| add | **1.9×** | 1.5e-32 | full DD |
| sub | **2.1×** | 6.1e-33 | full DD |
| mul | **8.6×** | 3.5e-32 | full DD |
| div | **3.1×** | 5.8e-32 | full DD |
| sqrt | **39×** | 5.4e-32 | full DD |
| add (mf+dp) | **2.6×** | exact | full DD |
| mul (dp\*mf) | **8.9×** | 2.7e-32 | full DD |

### Unary

| op | speedup | max\_rel | precision |
|---|---|---|---|
| abs | 1.4× | exact | exact |
| neg | **2.2×** | exact | exact |
| aint | 1.6× | exact | exact |
| anint | **2.4×** | exact | exact |
| fraction | 1.4× | exact | exact |
| scale | 1.0× | exact | exact |
| set\_exponent | 1.8× | exact | exact |

### Binary

| op | speedup | max\_rel | precision |
|---|---|---|---|
| min | **2.2×** | 6.1e-33 | full DD |
| max | **2.1×** | 5.9e-33 | full DD |
| min3 | **2.2×** | 5.9e-33 | full DD |
| max3 | **4.3×** | 6.0e-33 | full DD |
| sign | 1.2× | exact | exact |
| dim | **2.8×** | 6.1e-33 | full DD |
| hypot | **4.7×** | 8.1e-32 | full DD |
| mod | 0.53× | 1.6e-32 | full DD |
| modulo | **1.3×** | 1.6e-32 | full DD |

### Exponential / logarithmic

| op | speedup | max\_rel | precision |
|---|---|---|---|
| exp | **2.9×** | 4.9e-30 | full DD |
| log | **5.0×** | 3.1e-32 | full DD |
| log10 | **6.6×** | 3.1e-32 | full DD |
| pow | **4.3×** | 1.2e-30 | full DD |
| pow\_int | **16×** | 2.2e-32 | full DD |

### Trigonometric

| op | speedup | max\_rel | precision |
|---|---|---|---|
| sin | **2.7×** | 7.0e-25 | deriv-corrected |
| cos | **2.7×** | 9.9e-25 | deriv-corrected |
| tan | 1.3× | 9.9e-25 | deriv-corrected |
| asin | **2.0×** | 3.5e-32 | full DD |
| acos | 1.9× | 2.3e-32 | full DD |
| atan | 1.2× | 7.5e-23 | deriv-corrected |
| atan2 | 1.2× | 3.6e-32 | full DD |

### Hyperbolic

| op | speedup | max\_rel | precision |
|---|---|---|---|
| sinh | **2.3×** | 4.9e-30 | full DD |
| cosh | 1.6× | 4.9e-30 | full DD |
| tanh | **2.7×** | 1.5e-18 | deriv-corrected |
| asinh | **6.1×** | 2.4e-30 | full DD |
| acosh | **5.6×** | 3.5e-32 | full DD |
| atanh | **5.1×** | 1.1e-30 | full DD |

### Error / special functions

| op | speedup | max\_rel | precision |
|---|---|---|---|
| erf | **4.8×** | 3.6e-19 | deriv-corrected |
| erfc | **4.5×** | 4.5e-16 | deriv-corrected |
| erfc\_scaled | **162×** | 3.7e-16 | single-double |
| gamma | **112×** | 1.3e-16 | single-double |
| log\_gamma | **68×** | 3.7e-16 | single-double |
| bessel\_j0 | **96×** | 3.2e-16 | single-double |
| bessel\_j1 | **101×** | 2.5e-15 | single-double |
| bessel\_jn(3,.) | **104×** | 1.2e-15 | single-double |
| bessel\_y0 | **144×** | 1.2e-15 | single-double |
| bessel\_y1 | **145×** | 7.5e-16 | single-double |
| bessel\_yn(3,.) | **165×** | 5.2e-14 | single-double |

### Complex arithmetic

| op | speedup | max\_rel | precision |
|---|---|---|---|
| cx\_add | **3.4×** | 1.3e-32 | full DD |
| cx\_sub | **3.5×** | 5.9e-33 | full DD |
| cx\_mul | **7.4×** | 1.5e-32 | full DD |
| cx\_div | **7.3×** | 5.4e-32 (re) / 9.5e-17 (im) | full DD / deriv-corrected |
| cx\_conjg | 1.8× | exact | exact |
| cx\_abs | **4.0×** | 5.4e-32 | full DD |

### Complex transcendentals

| op | speedup | max\_rel | precision |
|---|---|---|---|
| cx\_sqrt | **5.8×** | 5.9e-32 | full DD |
| cx\_exp | **2.1×** | 1.5e-28 | deriv-corrected |
| cx\_log | **2.1×** | 1.0e-31 | full DD |
| cx\_sin | 1.9× | 3.9e-29 | deriv-corrected |
| cx\_cos | 1.9× | 5.8e-30 | full DD |
| cx\_tan | 0.96× | 9.4e-31 | full DD |
| cx\_sinh | **2.1×** | 1.5e-28 | deriv-corrected |
| cx\_cosh | **2.0×** | 3.7e-30 | full DD |
| cx\_tanh | 1.1× | 1.5e-30 | full DD |
| cx\_asin | **2.5×** | 4.3e-25 (re) / 6.2e-31 (im) | deriv-corrected |
| cx\_acos | **2.5×** | 1.2e-32 (re) / 6.2e-31 (im) | full DD |
| cx\_atan | 1.8× | 5.1e-32 (re) / 1.7e-31 (im) | full DD |
| cx\_asinh | **2.5×** | 5.4e-23 (re) / 7.2e-32 (im) | deriv-corrected |
| cx\_acosh | **2.2×** | 1.3e-30 (re) / 1.3e-32 (im) | full DD |
| cx\_atanh | **2.1×** | 2.5e-23 (re) / 4.3e-32 (im) | deriv-corrected |

### Array reductions

| op | speedup | max\_rel | precision |
|---|---|---|---|
| arr\_sum (n=8) | 1.2× | 6.5e-32 | full DD |
| arr\_product (n=8) | **3.4×** | 7.7e-53 | full DD |
| arr\_maxval (n=8) | **4.6×** | 5.7e-33 | full DD |
| arr\_minval (n=8) | **4.5×** | 5.9e-33 | full DD |
| arr\_dot (n=8) | **5.0×** | 2.5e-31 | full DD (fused FMA) |
| arr\_norm2 (n=8) | **5.6×** | 5.1e-32 | full DD |
| arr\_matmul (8×8\*8) | **2.1×** | 9.9e-31 | full DD (fused FMA) |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.
Precision characteristics are the same as the Fortran version (same
algorithms), so only speedup is shown.

### Arithmetic

| op | speedup |
|---|---|
| add | **3.0×** |
| sub | **4.4×** |
| mul | **11×** |
| div | **3.9×** |
| sqrt | **61×** |
| cbrt | **43×** |
| fma | **98×** |
| abs | **2.1×** |
| neg | **2.1×** |

### Rounding

| op | speedup |
|---|---|
| floor | **2.9×** |
| ceil | **3.2×** |
| trunc | **2.6×** |
| round | 1.0× |
| rint | **11×** |
| nearbyint | **34×** |

### Binary

| op | speedup |
|---|---|
| fmin | **5.2×** |
| fmax | **5.4×** |
| fdim | **6.5×** |
| copysign | 1.9× |
| fmod | 0.66× |
| hypot | **42×** |
| ldexp(.,5) | **2.1×** |

### Exponential / logarithmic

| op | speedup |
|---|---|
| exp | **3.1×** |
| exp2 | **3.6×** |
| expm1 | **4.6×** |
| log | **5.0×** |
| log10 | **6.7×** |
| log2 | **6.2×** |
| log1p | **5.5×** |
| pow | **4.5×** |

### Trigonometric

| op | speedup |
|---|---|
| sin | **3.2×** |
| cos | **3.1×** |
| tan | 1.4× |
| asin | **2.1×** |
| acos | **2.1×** |
| atan | 1.3× |
| atan2 | 1.1× |

### Hyperbolic

| op | speedup |
|---|---|
| sinh | **2.6×** |
| cosh | 1.8× |
| tanh | **2.8×** |
| asinh | **6.9×** |
| acosh | **6.3×** |
| atanh | **5.4×** |

### Error / special functions

| op | speedup |
|---|---|
| erf | **5.0×** |
| erfc | **5.2×** |
| tgamma | **105×** |
| lgamma | **71×** |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster. The DD `mod` uses a floor-multiple reduction loop
  (for small quotients) or a full DD divide chain (for large quotients);
  libquadmath's `fmodq` uses a specialized bit-level remainder algorithm.
  `modulo` (Fortran) now beats qp at 1.3× thanks to the iterative approach.
  Precision degrades as `~10^(log10(quotient) - 31)` for large quotients,
  which is the inherent DD precision limit.

- **Single-double precision** functions (gamma, bessel, erfc\_scaled, etc.)
  achieve 60–165× speedup by evaluating `f(hi)` via the leading-limb libm
  intrinsic without a lo-limb correction. This gives `double`-level accuracy
  (~15 digits) rather than full DD (~31 digits). For applications needing
  full DD precision on these functions, a polynomial or series expansion
  would be required (at significant implementation cost and reduced speedup).

- **Deriv-corrected** functions (sin, cos, tan, atan, erf, tanh, etc.)
  use `f(hi) + f'(hi) * lo` to recover most of the DD precision. The
  correction is exact to first order, giving ~20–25 digits in typical cases.

- The **Fortran** multifloats module uses `elemental` functions on a
  `sequence` derived type. gfortran's ABI passes/returns these via hidden
  pointers (not in FP registers), adding ~1.5× overhead vs the C++ header-
  only version even with LTO inlining. See the performance note at the top
  of `fsrc/multifloats.fypp` for details and the `bind(c)` escape hatch.

- **Array reductions** (dot\_product, matmul) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`, with periodic
  renormalization (configurable via `mf_set_fma_renorm_interval`).
