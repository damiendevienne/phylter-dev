#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace phylter2 {

class Matrix {
 public:
  Matrix() = default;
  Matrix(std::size_t rows, std::size_t cols, double value = 0.0);

  [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
  [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
  [[nodiscard]] bool square() const noexcept { return rows_ == cols_; }

  [[nodiscard]] double& operator()(std::size_t row, std::size_t col);
  [[nodiscard]] double operator()(std::size_t row, std::size_t col) const;
  [[nodiscard]] std::span<double> data() noexcept { return values_; }
  [[nodiscard]] std::span<const double> data() const noexcept { return values_; }

  [[nodiscard]] static Matrix multiply(const Matrix& left, const Matrix& right);
  [[nodiscard]] static Matrix identity(std::size_t n);

 private:
  std::size_t rows_{0};
  std::size_t cols_{0};
  std::vector<double> values_;
};

[[nodiscard]] double max_abs_difference(const Matrix& left, const Matrix& right);

}  // namespace phylter2
