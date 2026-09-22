#include "phylter2/lapack.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

extern "C" {
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, int* info);
}

namespace phylter2 {

SymmetricEigenResult symmetric_eigen(const Matrix& input) {
  if (!input.square() || input.rows() == 0) {
    throw std::invalid_argument("eigen decomposition requires a non-empty square matrix");
  }
  if (input.rows() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("matrix is too large for the LAPACK interface");
  }
  for (double value : input.data())
    if (!std::isfinite(value)) throw std::invalid_argument("eigen decomposition requires finite values");

  int n = static_cast<int>(input.rows());
  int lda = n;
  int info = 0;
  int workspace_size = -1;
  double workspace_query = 0.0;
  std::vector<double> lapack_matrix(input.rows() * input.cols());
  // Input is symmetric, so its row-major and column-major physical layouts agree.
  std::copy(input.data().begin(), input.data().end(), lapack_matrix.begin());
  std::vector<double> values(input.rows());
  char jobz = 'V';
  char uplo = 'U';
  dsyev_(&jobz, &uplo, &n, lapack_matrix.data(), &lda,
         values.data(), &workspace_query, &workspace_size, &info);
  if (info != 0 || !std::isfinite(workspace_query) || workspace_query < 1.0 ||
      workspace_query > std::numeric_limits<int>::max()) {
    throw std::runtime_error("LAPACK workspace query failed");
  }
  workspace_size = static_cast<int>(workspace_query);
  std::vector<double> workspace(static_cast<std::size_t>(workspace_size));
  dsyev_(&jobz, &uplo, &n, lapack_matrix.data(), &lda,
         values.data(), workspace.data(), &workspace_size, &info);
  if (info != 0) {
    throw std::runtime_error("LAPACK symmetric eigen decomposition failed");
  }

  std::vector<std::size_t> order(input.rows());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&values](const std::size_t left, const std::size_t right) {
    return values[left] > values[right];
  });
  SymmetricEigenResult result{.values = std::vector<double>(input.rows()),
                               .vectors = Matrix(input.rows(), input.cols())};
  for (std::size_t new_col = 0; new_col < order.size(); ++new_col) {
    const std::size_t old_col = order[new_col];
    result.values[new_col] = values[old_col];
    for (std::size_t row = 0; row < input.rows(); ++row) {
      result.vectors(row, new_col) = lapack_matrix[old_col * input.rows() + row];
    }
  }
  return result;
}

}  // namespace phylter2
