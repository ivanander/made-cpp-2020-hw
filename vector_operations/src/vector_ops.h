#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

namespace task {
constexpr double eps = 1e-7;
using VectorDouble = std::vector<double>;
using VectorInt = std::vector<int>;

VectorDouble operator+(const VectorDouble &lhs, const VectorDouble &rhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(rhs),
                 begin(res),
                 std::plus<double>{});
  return res;
}

VectorDouble operator+(const VectorDouble &lhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(res),
                 [](double val) { return +val; });
  return res;
}

VectorDouble operator-(const VectorDouble &lhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(res),
                 [](double val) { return -val; });
  return res;
}

VectorDouble operator-(const VectorDouble &lhs, const VectorDouble &rhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(rhs),
                 begin(res),
                 std::minus<double>{});
  return res;
}

double operator*(const VectorDouble &lhs, const VectorDouble &rhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(rhs),
                 begin(res),
                 std::multiplies<double>{});
  return std::accumulate(begin(res), end(res), 0.0);
}

double length(const VectorDouble &vec) {
  double res = 0;
  for (double val: vec) {
    res += val * val;
  }
  return sqrt(res);
}

bool operator||(const VectorDouble &lhs, const VectorDouble &rhs) {
  double scalar_mul = lhs * rhs;
  double length_mul = length(lhs) * length(rhs);
  return abs(abs(scalar_mul) - length_mul) < eps;
}

bool operator&&(const VectorDouble &lhs, const VectorDouble &rhs) {
  double scalar_mul = lhs * rhs;
  double length_mul = length(lhs) * length(rhs);
  return abs(scalar_mul - length_mul) < eps;
}

// only for 3-dimensional vectors
VectorDouble operator%(const VectorDouble &lhs, const VectorDouble &rhs) {
  auto res = lhs;
  res[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
  res[1] = -(lhs[0] * rhs[2] - lhs[2] * rhs[0]);
  res[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
  return res;
}

std::istream &operator>>(std::istream &input, VectorDouble &lhs) {
  size_t size;
  input >> size;
  lhs.resize(size);
  for (size_t i = 0; i < size; ++i)
    input >> lhs[i];
  return input;
}

std::ostream &operator<<(std::ostream &output, const VectorDouble &lhs) {
  for (int i = 0; i < lhs.size() - 1; ++i)
    output << lhs[i] << " ";
  if (!lhs.empty())
    output << lhs.back() << "\n";
  return output;
}

void reverse(VectorDouble &lhs) {
  auto left = lhs.begin();
  auto right = lhs.end() - 1;
  for (; left < right; ++left, --right) {
    std::swap(*left, *right);
  }
}

VectorInt operator|(const VectorInt &lhs, const VectorInt &rhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(rhs),
                 begin(res),
                 [](int lhs, int rhs) { return lhs | rhs; }
  );
  return res;
}

VectorInt operator&(const VectorInt &lhs, const VectorInt &rhs) {
  auto res = lhs;
  std::transform(begin(lhs),
                 end(lhs),
                 begin(rhs),
                 begin(res),
                 [](int lhs, int rhs) { return lhs & rhs; });
  return res;
}
}  // namespace task
