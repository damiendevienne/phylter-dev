#include "phylter2/rv.hpp"
#include <cmath>
#include <stdexcept>
namespace phylter2 {
namespace {
double inner(const Matrix& left, const Matrix& right) {
  const auto a = left.data(), b = right.data();
  const auto n = left.rows();
  double result = 0;
  for (std::size_t row = 0; row < n; ++row) {
    result += a[row*n+row]*b[row*n+row];
    for (std::size_t col = row+1; col < n; ++col)
      result += 2*a[row*n+col]*b[row*n+col];
  }
  return result;
}
}
RvOperator::RvOperator(const std::vector<Matrix>& centered) : centered_(centered) {
  if (centered.empty()) throw std::invalid_argument("empty RV input");
  const auto n = centered.front().rows();
  if (!n) throw std::invalid_argument("empty RV matrix");
  for (const auto& matrix : centered) {
    if (!matrix.square() || matrix.rows()!=n) throw std::invalid_argument("incompatible RV matrices");
    const double length = std::sqrt(inner(matrix,matrix));
    if (!std::isfinite(length) || length==0) throw std::runtime_error("a centered matrix has zero or non-finite norm");
    norms_.push_back(length);
  }
}
std::vector<double> RvOperator::operator()(const std::vector<double>& x) const {
  const auto genes = centered_.size(), n = centered_.front().rows();
  if (x.size()!=genes) throw std::invalid_argument("incompatible RV vector");
  Matrix aggregate(n,n);
  auto buffer = aggregate.data();
  std::vector<double> scaled(genes);
  for (std::size_t gene = 0; gene < genes; ++gene) scaled[gene] = x[gene]/norms_[gene];
  // Each row is owned by one thread. Gene accumulation order stays fixed, so
  // different thread counts do not change reduction order or convergence.
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(n); ++r) {
    const auto row = static_cast<std::size_t>(r);
    for (std::size_t gene = 0; gene < genes; ++gene) {
      const auto matrix = centered_[gene].data();
      for (std::size_t col = row; col < n; ++col)
        buffer[row*n+col] += scaled[gene]*matrix[row*n+col];
    }
  }
  std::vector<double> y(genes);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (std::ptrdiff_t g = 0; g < static_cast<std::ptrdiff_t>(genes); ++g) {
    const auto gene = static_cast<std::size_t>(g);
    y[gene] = inner(centered_[gene],aggregate)/norms_[gene];
  }
  return y;
}
Matrix RvOperator::materialize() const {
  const auto genes = centered_.size();
  Matrix rv(genes,genes);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(genes); ++r) {
    const auto row = static_cast<std::size_t>(r);
    for (std::size_t col = row; col < genes; ++col)
      rv(row,col) = rv(col,row) = inner(centered_[row],centered_[col])/(norms_[row]*norms_[col]);
  }
  return rv;
}
}
