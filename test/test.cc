#include "matmul.hh"
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <typeinfo>
#include <vector>

template <typename T>
void test_matmul_square(std::size_t n, std::size_t seed = 0) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> dist;

  auto generator = [&rng, &dist]() -> T { return dist(rng); };

  std::vector<T> a;
  a.reserve(n * n);
  std::generate_n(std::back_inserter(a), n, generator);

  std::vector<T> b;
  b.reserve(n * n);
  std::generate_n(std::back_inserter(b), n, generator);

  std::vector<T> c(n * n);

  auto start_time = std::chrono::system_clock::now();
  matmul(n, n, n, a.data(), n, b.data(), n, c.data(), n);
  auto end_time = std::chrono::system_clock::now();
  // TODO: do something about the duration
  std::cout << typeid(T).name() << "\t" << n << "\t" << seed << "\t"
            << std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                     start_time)
                       .count() *
                   1e-3
            << "ms\n";
}

int main() { test_matmul_square<double>(512, 0); }
