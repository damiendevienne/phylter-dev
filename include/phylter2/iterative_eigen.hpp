#pragma once
#include <cstddef>
#include <functional>
#include <vector>
namespace phylter2 {
using SymmetricOperator = std::function<std::vector<double>(const std::vector<double>&)>;
struct KrylovOptions {
  // A 1e-12 residual is stable across LAPACK implementations. Requiring a
  // value near machine precision can make an otherwise converged solve fail
  // on platforms with different rounding behavior (notably Apple Accelerate).
  double tolerance{1e-12};
  std::size_t max_basis{48};
  std::size_t max_products{2048};
};
struct LeadingEigenpair {
  double value{};
  std::vector<double> vector;
  double relative_residual{};
  std::size_t products{};
};
// Largest algebraic eigenpair of a symmetric operator; for PSD Gram matrices
// this is also the largest-magnitude eigenpair. Bounded, fully reorthogonalized
// Krylov space with explicit Ritz-vector restarts. Throws on non-convergence.
LeadingEigenpair leading_eigenpair(std::size_t size, const SymmetricOperator& apply,
                                  const KrylovOptions& options = {});
}
