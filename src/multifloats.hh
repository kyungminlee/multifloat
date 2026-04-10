#pragma once

#include <array>

template <typename T, std::size_t N> struct MultiFloat {
  static_assert(N > 0, "N must be positive");
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
      if ((_limbs[i] != rhs._limbs[i])) {
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

  constexpr MultiFloat &operator+=(MultiFloat const &rhs);

  constexpr MultiFloat &operator-=(MultiFloat const &rhs);

  constexpr MultiFloat operator+(MultiFloat const &rhs) const;

  constexpr MultiFloat operator-(MultiFloat const &rhs) const;
};

template <typename T>
constexpr void two_sum(T const &a, T const &b, T &sum, T &err) {
  T s = a + b;
  T a_prime = s - b;
  T b_prime = s - a_prime;
  T a_err = a - a_prime;
  T b_err = b - b_prime;
  err = a_err + b_err;
  sum = s;
}

template <typename T>
constexpr void fast_two_sum(T const &a, T const &b, T &sum, T &err) {
  T s = a + b;
  T b_prime = s - a;
  err = b - b_prime;
  sum = s;
}

template <typename T>
constexpr void mfadd(T const (&x)[1], T const (&y)[1], T const (&z)[1]) {
  z[0] = x[0] + y[0];
}

template <typename T>
constexpr void mfadd(T const (&x)[2], T const (&y)[2], T const (&z)[2]) {
  T s, e;
  two_sum(x[0], y[0], s, e);
  e += x[1] + y[1];
  fast_two_sum(s, e, z[0], z[1]);
}

// TODO implement fma

template <typename T>
constexpr void two_prod(T const &a, T const &b, T &prod, T &err) {
  T p = a * b;
  T e = fma(a, b, -p);
  prod = p;
  err = e;
}

template <typename T>
constexpr void mfmul(T const (&x)[1], T const (&y)[1], T const (&z)[1]) {
  z[0] = x[0] * y[0];
}

template <typename T>
constexpr void mfmul(T const (&x)[2], T const (&y)[2], T const (&z)[2]) {
  T p, e;
  two_prod(x[0], y[0], p, e);
  e += x[0] * y[1];
  e += x[1] * y[0];
  e += x[1] * y[1];
  fast_two_sum(p, e, z[0], z[1]);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> &
MultiFloat<T, N>::operator+=(MultiFloat<T, N> const &rhs) {
  for (std::size_t i = 0; i < N; ++i) {
    mfadd(_limbs, rhs._limbs, _limbs);
  }
  return *this;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> &
MultiFloat<T, N>::operator-=(MultiFloat<T, N> const &rhs) {
  return (*this) += (-rhs);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N>
MultiFloat<T, N>::operator+(MultiFloat<T, N> const &rhs) const {
  MultiFloat<T, N> out(*this);
  out += rhs;
  return out;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N>
MultiFloat<T, N>::operator-(MultiFloat<T, N> const &rhs) const {
  MultiFloat<T, N> out(*this);
  out -= rhs;
  return out;
}

// min
// max
// signbit
// exponent
// issubnormal
// isfinite
// isinf
// isone
// iszero
// isinteger
// ldexp
//
