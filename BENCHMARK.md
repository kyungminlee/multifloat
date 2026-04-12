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

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. The "speedup" column is `qp_time / mf_time`:
values > 1× mean multifloats is faster.

### Arithmetic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| add | 0.0033 | 0.0013 | **2.4×** |
| sub | 0.0031 | 0.0014 | **2.2×** |
| mul | 0.0069 | 0.0007 | **9.3×** |
| div | 0.0109 | 0.0036 | **3.0×** |
| sqrt | 0.1538 | 0.0039 | **39×** |
| add (mf+dp) | 0.0027 | 0.0011 | **2.4×** |
| mul (dp\*mf) | 0.0071 | 0.0008 | **9.4×** |

### Unary

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| abs | 0.0012 | 0.0008 | 1.5× |
| neg | 0.0012 | 0.0005 | **2.2×** |
| aint | 0.0020 | 0.0012 | 1.7× |
| anint | 0.0022 | 0.0010 | **2.1×** |
| fraction | 0.0032 | 0.0023 | 1.4× |
| scale | 0.0018 | 0.0018 | 1.0× |
| set\_exponent | 0.0037 | 0.0022 | 1.7× |

### Binary

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| min | 0.0024 | 0.0010 | **2.3×** |
| max | 0.0023 | 0.0010 | **2.4×** |
| min3 | 0.0047 | 0.0014 | **3.3×** |
| max3 | 0.0031 | 0.0014 | **2.2×** |
| sign | 0.0015 | 0.0013 | 1.2× |
| dim | 0.0066 | 0.0030 | **2.2×** |
| hypot | 0.2041 | 0.0441 | **4.6×** |
| mod | 0.0026 | 0.0180 | 0.15× |
| modulo | 0.0081 | 0.0225 | 0.36× |

### Exponential / logarithmic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| exp | 0.0221 | 0.0078 | **2.9×** |
| log | 0.0275 | 0.0055 | **5.0×** |
| log10 | 0.0364 | 0.0056 | **6.5×** |
| pow | 0.0694 | 0.0167 | **4.1×** |
| pow\_int | 0.0023 | 0.0001 | **16×** |

### Trigonometric

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sin | 0.0222 | 0.0081 | **2.7×** |
| cos | 0.0227 | 0.0081 | **2.8×** |
| tan | 0.0208 | 0.0168 | 1.2× |
| asin | 0.0375 | 0.0187 | **2.0×** |
| acos | 0.0379 | 0.0204 | 1.9× |
| atan | 0.0233 | 0.0197 | 1.2× |
| atan2 | 0.0261 | 0.0219 | 1.2× |

### Hyperbolic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sinh | 0.0374 | 0.0160 | **2.3×** |
| cosh | 0.0258 | 0.0167 | 1.5× |
| tanh | 0.0367 | 0.0132 | **2.8×** |
| asinh | 0.0543 | 0.0089 | **6.1×** |
| acosh | 0.0488 | 0.0088 | **5.5×** |
| atanh | 0.0434 | 0.0082 | **5.3×** |

### Error / special functions

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| erf | 0.0071 | 0.0013 | **5.4×** |
| erfc | 0.0070 | 0.0014 | **5.0×** |
| erfc\_scaled | 0.0094 | 0.0001 | **163×** |
| gamma | 0.0111 | 0.0001 | **114×** |
| log\_gamma | 0.0033 | 0.0000 | **67×** |
| bessel\_j0 | 0.0095 | 0.0001 | **94×** |
| bessel\_j1 | 0.0093 | 0.0001 | **96×** |
| bessel\_jn(3,.) | 0.0284 | 0.0003 | **102×** |
| bessel\_y0 | 0.0143 | 0.0001 | **142×** |
| bessel\_y1 | 0.0144 | 0.0001 | **145×** |
| bessel\_yn(3,.) | 0.0284 | 0.0002 | **163×** |

### Complex arithmetic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| cx\_add | 0.0087 | 0.0027 | **3.3×** |
| cx\_sub | 0.0095 | 0.0026 | **3.7×** |
| cx\_mul | 0.0288 | 0.0041 | **7.1×** |
| cx\_div | 0.0940 | 0.0129 | **7.3×** |
| cx\_conjg | 0.0020 | 0.0011 | 1.9× |
| cx\_abs | 0.2007 | 0.0500 | **4.0×** |

### Complex transcendentals

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| cx\_sqrt | 0.0042 | 0.0008 | **5.6×** |
| cx\_exp | 0.0052 | 0.0025 | **2.1×** |
| cx\_log | 0.0082 | 0.0036 | **2.3×** |
| cx\_sin | 0.0092 | 0.0051 | 1.8× |
| cx\_cos | 0.0089 | 0.0050 | 1.8× |
| cx\_tan | 0.0098 | 0.0104 | 0.94× |
| cx\_sinh | 0.0099 | 0.0048 | **2.0×** |
| cx\_cosh | 0.0096 | 0.0048 | **2.0×** |
| cx\_tanh | 0.0102 | 0.0099 | 1.0× |
| cx\_asin | 0.0123 | 0.0048 | **2.6×** |
| cx\_acos | 0.0123 | 0.0048 | **2.6×** |
| cx\_atan | 0.0071 | 0.0040 | 1.8× |
| cx\_asinh | 0.0122 | 0.0049 | **2.5×** |
| cx\_acosh | 0.0119 | 0.0056 | **2.1×** |
| cx\_atanh | 0.0083 | 0.0039 | **2.1×** |

### Array reductions

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| arr\_sum (n=8) | 0.0275 | 0.0234 | 1.2× |
| arr\_product (n=8) | 0.0409 | 0.0117 | **3.5×** |
| arr\_maxval (n=8) | 0.0089 | 0.0019 | **4.6×** |
| arr\_minval (n=8) | 0.0089 | 0.0019 | **4.8×** |
| arr\_dot (n=8) | 0.0337 | 0.0061 | **5.5×** |
| arr\_norm2 (n=8) | 0.1577 | 0.0282 | **5.6×** |
| arr\_matmul (8×8\*8) | 0.0299 | 0.0151 | **2.0×** |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.

### Arithmetic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| add | 0.0037 | 0.0010 | **3.7×** |
| sub | 0.0051 | 0.0011 | **4.5×** |
| mul | 0.0082 | 0.0008 | **10.8×** |
| div | 0.0128 | 0.0032 | **4.0×** |
| sqrt | 0.1538 | 0.0026 | **58×** |
| cbrt | 0.2319 | 0.0053 | **44×** |
| fma | 0.1015 | 0.0011 | **95×** |
| abs | 0.0015 | 0.0007 | **2.1×** |
| neg | 0.0011 | 0.0005 | **2.2×** |

### Rounding

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| floor | 0.0020 | 0.0006 | **3.1×** |
| ceil | 0.0021 | 0.0007 | **3.1×** |
| trunc | 0.0021 | 0.0008 | **2.6×** |
| round | 0.0021 | 0.0021 | 1.0× |
| rint | 0.0067 | 0.0007 | **10×** |
| nearbyint | 0.0228 | 0.0007 | **35×** |

### Binary

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| fmin | 0.0048 | 0.0009 | **5.4×** |
| fmax | 0.0046 | 0.0008 | **5.6×** |
| fdim | 0.0056 | 0.0009 | **6.0×** |
| copysign | 0.0015 | 0.0008 | 1.9× |
| fmod | 0.0028 | 0.0081 | 0.34× |
| hypot | 0.2103 | 0.0054 | **39×** |
| ldexp(.,5) | 0.0039 | 0.0018 | **2.2×** |

### Exponential / logarithmic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| exp | 0.0258 | 0.0071 | **3.6×** |
| exp2 | 0.0255 | 0.0070 | **3.7×** |
| expm1 | 0.0343 | 0.0074 | **4.6×** |
| log | 0.0289 | 0.0058 | **5.0×** |
| log10 | 0.0373 | 0.0057 | **6.6×** |
| log2 | 0.0441 | 0.0057 | **7.8×** |
| log1p | 0.0341 | 0.0058 | **5.8×** |
| pow | 0.0745 | 0.0388 | 1.9× |

### Trigonometric

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sin | 0.0394 | 0.0075 | **5.3×** |
| cos | 0.0219 | 0.0073 | **3.0×** |
| tan | 0.0205 | 0.0149 | 1.4× |
| asin | 0.0364 | 0.0178 | **2.1×** |
| acos | 0.0378 | 0.0180 | **2.1×** |
| atan | 0.0228 | 0.0179 | 1.3× |
| atan2 | 0.0256 | 0.0227 | 1.1× |

### Hyperbolic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sinh | 0.0376 | 0.0143 | **2.6×** |
| cosh | 0.0252 | 0.0145 | 1.7× |
| tanh | 0.0378 | 0.0132 | **2.9×** |
| asinh | 0.0535 | 0.0077 | **7.0×** |
| acosh | 0.0479 | 0.0078 | **6.2×** |
| atanh | 0.0426 | 0.0078 | **5.4×** |

### Error / special functions

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| erf | 0.0073 | 0.0013 | **5.4×** |
| erfc | 0.0070 | 0.0014 | **5.0×** |
| tgamma | 0.0108 | 0.0001 | **114×** |
| lgamma | 0.0032 | 0.0000 | **66×** |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster (~3–4×). The DD `mod` requires a full DD divide +
  truncate + multiply + subtract chain; libquadmath's `fmodq` uses a
  specialized bit-level remainder algorithm with no DD equivalent.

- The **Fortran** multifloats module uses `elemental` functions on a
  `sequence` derived type. gfortran's ABI passes/returns these via hidden
  pointers (not in FP registers), adding ~1.5× overhead vs the C++ header-
  only version even with LTO inlining. See the performance note at the top
  of `fsrc/multifloats.fypp` for details and the `bind(c)` escape hatch.

- **Special functions** (gamma, bessel, erfc\_scaled) achieve 60–163×
  speedup because the DD kernels use a leading-limb libm call + derivative
  correction, which is inherently cheaper than libquadmath's full software
  quad implementation of these functions.

- **Array reductions** (dot\_product, matmul) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`, with periodic
  renormalization (configurable via `mf_set_fma_renorm_interval`).
