#include "phylter2/matrix.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

extern "C" void dgemm_(const char*, const char*, const int*, const int*, const int*,
                       const double*, const double*, const int*, const double*, const int*,
                       const double*, double*, const int*);

namespace phylter2 {

Matrix::Matrix(const std::size_t rows, const std::size_t cols, const double value)
    : rows_(rows), cols_(cols) {
  if (cols && rows > values_.max_size() / cols) throw std::length_error("matrix dimensions overflow");
  values_.assign(rows * cols, value);
}

double& Matrix::operator()(const std::size_t row, const std::size_t col) {
  if (row >= rows_ || col >= cols_) {
    throw std::out_of_range("matrix index out of range");
  }
  return values_[row * cols_ + col];
}

double Matrix::operator()(const std::size_t row, const std::size_t col) const {
  if (row >= rows_ || col >= cols_) {
    throw std::out_of_range("matrix index out of range");
  }
  return values_[row * cols_ + col];
}

Matrix Matrix::multiply(const Matrix& left, const Matrix& right) {
  if (left.cols_ != right.rows_) {
    throw std::invalid_argument("incompatible matrix dimensions");
  }
  Matrix result(left.rows_, right.cols_);
  if (result.data().empty() || left.cols_ == 0) return result;
  const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (left.rows_ > limit || left.cols_ > limit || right.cols_ > limit)
    throw std::length_error("matrix too large for LP64 BLAS");
  const int m = static_cast<int>(right.cols_), n = static_cast<int>(left.rows_),
            k = static_cast<int>(left.cols_);
  const char normal = 'N';
  const double one = 1, zero = 0;
  // (AB)^T = B^T A^T allows direct use of contiguous row-major storage.
  dgemm_(&normal, &normal, &m, &n, &k, &one, right.values_.data(), &m,
         left.values_.data(), &k, &zero, result.values_.data(), &m);
  return result;
}

Matrix Matrix::identity(const std::size_t n) {
  Matrix result(n, n);
  for (std::size_t i = 0; i < n; ++i) {
    result(i, i) = 1.0;
  }
  return result;
}

double max_abs_difference(const Matrix& left, const Matrix& right) {
  if (left.rows() != right.rows() || left.cols() != right.cols()) {
    throw std::invalid_argument("incompatible matrix dimensions");
  }
  double result = 0.0;
  for (std::size_t i = 0; i < left.data().size(); ++i) {
    if (!std::isfinite(left.data()[i]) || !std::isfinite(right.data()[i]))
      return std::numeric_limits<double>::infinity();
    result = std::max(result, std::abs(left.data()[i] - right.data()[i]));
  }
  return result;
}

}  // namespace phylter2
