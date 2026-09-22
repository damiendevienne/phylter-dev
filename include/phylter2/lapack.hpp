#pragma once

#include "phylter2/matrix.hpp"

#include <vector>

namespace phylter2 {

struct SymmetricEigenResult {
  // Eigenvalues are in descending order; eigenvectors are stored by column.
  std::vector<double> values;
  Matrix vectors;
};

[[nodiscard]] SymmetricEigenResult symmetric_eigen(const Matrix& input);

}  // namespace phylter2
