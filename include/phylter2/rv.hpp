#pragma once
#include "phylter2/matrix.hpp"
#include "phylter2/iterative_eigen.hpp"
namespace phylter2 {
enum class RvMethod { matrix_free, dense };
struct RvOptions {
  RvMethod method{RvMethod::matrix_free};
  bool materialize{true}; // Direct DISTATIS API compatibility; CLI defaults false.
};
// References existing centered matrices; no additional gene-by-feature copy.
// Uses diagonal + twice upper triangle exactly as the R RV definition does,
// including when the legacy imputation fallback produces asymmetry.
class RvOperator {
 public:
  explicit RvOperator(const std::vector<Matrix>& centered);
  std::vector<double> operator()(const std::vector<double>& vector) const;
  Matrix materialize() const;
 private:
  const std::vector<Matrix>& centered_;
  std::vector<double> norms_;
};
}
