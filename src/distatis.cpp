#include "phylter2/distatis.hpp"

#include "phylter2/lapack.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace phylter2 {
namespace {

void validate_matrices(const std::vector<Matrix>& matrices) {
  if (matrices.size() < 2) {
    throw std::invalid_argument("DISTATIS requires at least two matrices");
  }
  const std::size_t species = matrices.front().rows();
  if (species < 3 || !matrices.front().square()) {
    throw std::invalid_argument("DISTATIS requires square matrices for at least three taxa");
  }
  for (const Matrix& matrix : matrices) {
    if (!matrix.square() || matrix.rows() != species) {
      throw std::invalid_argument("all DISTATIS matrices must have the same square dimension");
    }
    for (const double value : matrix.data()) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument("DISTATIS matrices must contain finite values");
      }
    }
  }
}

Matrix double_center(const Matrix& input) {
  const std::size_t n = input.rows();
  std::vector<double> row_mean(n, 0.0);
  std::vector<double> col_mean(n, 0.0);
  double overall = 0.0;
  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col) {
      row_mean[row] += input(row, col);
      col_mean[col] += input(row, col);
      overall += input(row, col);
    }
  }
  const double inverse_n = 1.0 / static_cast<double>(n);
  for (double& value : row_mean) value *= inverse_n;
  for (double& value : col_mean) value *= inverse_n;
  overall *= inverse_n * inverse_n;

  Matrix result(n, n);
  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col) {
      result(row, col) = -0.5 * (input(row, col) - col_mean[col] - row_mean[row] + overall);
    }
  }
  return result;
}

std::size_t choose_factors(const std::vector<double>& values, const std::size_t available) {
  double total = std::accumulate(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(available), 0.0);
  if (total == 0.0 || !std::isfinite(total)) {
    throw std::runtime_error("cannot choose factors from a zero or non-finite spectrum");
  }
  for (std::size_t index = 0; index < available; ++index) {
    double broken_stick = 0.0;
    for (std::size_t term = index + 1; term <= available; ++term) {
      broken_stick += 1.0 / static_cast<double>(term);
    }
    broken_stick /= static_cast<double>(available);
    if (values[index] / total <= broken_stick) {
      return std::max<std::size_t>(2, index);
    }
  }
  return available;
}

}  // namespace

DistatisResult distatis(const std::vector<Matrix>& matrices, const std::size_t requested_factors,
                        const double minimum_quality, const RvOptions& rv_options) {
  validate_matrices(matrices);
  const std::size_t genes = matrices.size();
  const std::size_t species = matrices.front().rows();
  std::vector<Matrix> centered;
  centered.reserve(genes);
  for (const Matrix& matrix : matrices) centered.push_back(double_center(matrix));

  const RvOperator rv_operator(centered);
  Matrix rv;
  LeadingEigenpair leading;
  if (rv_options.method == RvMethod::dense) {
    rv = rv_operator.materialize();
    const auto eig = symmetric_eigen(rv);
    leading.value = eig.values.front();
    leading.relative_residual = std::numeric_limits<double>::quiet_NaN(); // not evaluated
    leading.vector.resize(genes);
    for (std::size_t gene = 0; gene < genes; ++gene) leading.vector[gene] = eig.vectors(gene,0);
  } else {
    leading = leading_eigenpair(genes,[&](const auto& x) { return rv_operator(x); });
  }
  std::vector<double> alpha(genes);
  const double first_vector_sum = std::accumulate(leading.vector.begin(),leading.vector.end(),0.0);
  if (first_vector_sum == 0.0 || !std::isfinite(first_vector_sum))
    throw std::runtime_error("leading RV eigenvector cannot be normalized");
  for (std::size_t gene = 0; gene < genes; ++gene) alpha[gene] = leading.vector[gene]/first_vector_sum;
  const double quality = leading.value/static_cast<double>(genes);
  // Rejected proposals need neither compromise eigensystem nor projections.
  if (quality < minimum_quality) {
    DistatisResult rejected;
    rejected.quality = quality;
    return rejected;
  }
  // RV output is optional and is never needed to accept or reject a proposal.
  if (rv_options.materialize && rv.rows()==0) rv = rv_operator.materialize();
  if (!rv_options.materialize) rv = Matrix();
  Matrix compromise_centered(species, species);
  for (std::size_t gene = 0; gene < genes; ++gene) {
    for (std::size_t row = 0; row < species; ++row) {
      for (std::size_t col = 0; col < species; ++col) {
        compromise_centered(row, col) += alpha[gene] * centered[gene](row, col);
      }
    }
  }
  Matrix compromise(species, species);
  for (std::size_t row = 0; row < species; ++row) {
    for (std::size_t col = 0; col < species; ++col) {
      compromise(row, col) = -2.0 * compromise_centered(row, col) +
                             compromise_centered(row, row) + compromise_centered(col, col);
    }
  }

  const SymmetricEigenResult full_eigen = symmetric_eigen(compromise_centered);
  // RSpectra which="LM": select by magnitude, return in algebraic order.
  const std::size_t selected = requested_factors == 0 ? species - 1 : requested_factors;
  if (selected == 0 || selected >= species)
    throw std::invalid_argument("requested number of factors is outside the supported range");
  std::vector<std::size_t> axes(species);
  std::iota(axes.begin(), axes.end(), 0);
  std::stable_sort(axes.begin(), axes.end(), [&](auto a, auto b) {
    return std::abs(full_eigen.values[a]) > std::abs(full_eigen.values[b]);
  });
  axes.resize(selected);
  std::sort(axes.begin(), axes.end());
  SymmetricEigenResult compromise_eigen{.values = {}, .vectors = Matrix(species, selected)};
  for (std::size_t a = 0; a < selected; ++a) {
    compromise_eigen.values.push_back(full_eigen.values[axes[a]]);
    for (std::size_t i = 0; i < species; ++i)
      compromise_eigen.vectors(i,a) = full_eigen.vectors(i,axes[a]);
  }
  const std::size_t available = species - 1;
  const std::size_t factors = requested_factors == 0 ? choose_factors(compromise_eigen.values, available)
                                                     : requested_factors;
  if (factors == 0 || factors > available) {
    throw std::invalid_argument("requested number of factors is outside the supported range");
  }
  Matrix coordinates(species, factors);
  Matrix projection(species, factors);
  for (std::size_t axis = 0; axis < factors; ++axis) {
    const double singular = std::sqrt(std::abs(compromise_eigen.values[axis]));
    if (singular == 0.0 || !std::isfinite(singular)) {
      throw std::runtime_error("cannot project a zero or non-finite DISTATIS axis");
    }
    for (std::size_t species_index = 0; species_index < species; ++species_index) {
      coordinates(species_index, axis) = compromise_eigen.vectors(species_index, axis) * singular;
      projection(species_index, axis) = compromise_eigen.vectors(species_index, axis) / singular;
    }
  }
  std::vector<Matrix> partial;
  partial.reserve(genes);
  for (const Matrix& matrix : centered) partial.push_back(Matrix::multiply(matrix, projection));

  return {.coordinates = std::move(coordinates),
          .partial_coordinates = std::move(partial),
          .alpha = std::move(alpha),
          .eigenvalues = compromise_eigen.values,
          .rv = std::move(rv),
          .compromise = std::move(compromise),
          .quality = quality,
          .rv_products = leading.products,
          .rv_relative_residual = leading.relative_residual,
          .centered = std::move(centered)};
}

Matrix dist2wr(const DistatisResult& result) {
  const std::size_t species = result.coordinates.rows();
  const std::size_t genes = result.partial_coordinates.size();
  Matrix result_matrix(species, genes);
  for (std::size_t gene = 0; gene < genes; ++gene) {
    if (result.partial_coordinates[gene].rows() != species ||
        result.partial_coordinates[gene].cols() != result.coordinates.cols()) {
      throw std::invalid_argument("incompatible partial DISTATIS coordinates");
    }
    for (std::size_t species_index = 0; species_index < species; ++species_index) {
      double squared_distance = 0.0;
      for (std::size_t axis = 0; axis < result.coordinates.cols(); ++axis) {
        const double difference = result.partial_coordinates[gene](species_index, axis) -
                                  result.coordinates(species_index, axis);
        squared_distance += difference * difference;
      }
      result_matrix(species_index, gene) = std::sqrt(squared_distance);
    }
  }
  return result_matrix;
}

}  // namespace phylter2
