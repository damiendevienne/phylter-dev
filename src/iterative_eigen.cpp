#include "phylter2/iterative_eigen.hpp"
#include "phylter2/lapack.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
namespace phylter2 {
namespace {
using Vector = std::vector<double>;
using Basis = std::vector<Vector>;
double dot(const Vector& a, const Vector& b) {
  long double sum = 0;
  for (std::size_t i = 0; i < a.size(); ++i) sum += static_cast<long double>(a[i])*b[i];
  return static_cast<double>(sum);
}
double norm(const Vector& v) { return std::sqrt(dot(v,v)); }
bool append_orthogonal(Basis& basis, Vector v) {
  const double original = norm(v);
  if (!std::isfinite(original)) throw std::runtime_error("non-finite Krylov vector");
  // Two modified Gram-Schmidt passes suppress loss of orthogonality, including
  // rank-deficient RV operators with many repeated gene matrices.
  for (int pass = 0; pass < 2; ++pass)
    for (const auto& q : basis) {
      const double projection = dot(q,v);
      for (std::size_t i = 0; i < v.size(); ++i) v[i] -= projection*q[i];
    }
  const double length = norm(v);
  if (length <= 64*std::numeric_limits<double>::epsilon()*original || length == 0) return false;
  for (auto& x : v) x /= length;
  basis.push_back(std::move(v));
  return true;
}
Vector combination(const Basis& basis, const Matrix& coefficients, std::size_t column) {
  Vector v(basis.front().size());
  for (std::size_t j = 0; j < coefficients.rows(); ++j)
    for (std::size_t i = 0; i < v.size(); ++i) v[i] += coefficients(j,column)*basis[j][i];
  return v;
}
}
LeadingEigenpair leading_eigenpair(const std::size_t size, const SymmetricOperator& apply,
                                  const KrylovOptions& options) {
  if (!size || !std::isfinite(options.tolerance) || options.tolerance <= 0 ||
      options.max_basis < 4 || options.max_products < 1)
    throw std::invalid_argument("invalid iterative eigensolver options");
  const auto capacity = std::min(size,options.max_basis);
  Basis q, aq;
  std::size_t products = 0;
  auto checked_apply = [&](const Vector& v) {
    if (products == options.max_products)
      throw std::runtime_error("RV eigensolver did not converge; retry with --rv-method dense");
    ++products;
    auto result = apply(v);
    if (result.size() != size) throw std::runtime_error("operator returned an incompatible vector");
    for (double x : result)
      if (!std::isfinite(x)) throw std::runtime_error("operator returned a non-finite value");
    return result;
  };
  // A constant start is effective for correlated genes. A second, independent
  // deterministic pseudorandom start also probes eigenspaces orthogonal to it.
  std::uint64_t seed = 0x3c79ac492ba7b653ULL;
  auto random_vector = [&]() {
    Vector v(size);
    for (auto& value : v) {
      seed += 0x9e3779b97f4a7c15ULL;
      auto x = seed;
      x = (x ^ (x >> 30))*0xbf58476d1ce4e5b9ULL;
      x = (x ^ (x >> 27))*0x94d049bb133111ebULL;
      x ^= x >> 31;
      value = static_cast<double>(x >> 11)*0x1.0p-53 - .5;
    }
    return v;
  };
  append_orthogonal(q,Vector(size,1));
  if (size > 1) append_orthogonal(q,random_vector());
  while (true) {
    // Expand both starting directions, or the Ritz vectors retained on restart.
    while (aq.size() < q.size() && aq.size() < capacity) {
      auto image = checked_apply(q[aq.size()]);
      aq.push_back(image);
      if (q.size() < capacity) append_orthogonal(q,std::move(image));
      const auto m = aq.size();
      const bool exhausted = m == q.size();
      if (!exhausted && m < capacity && (m < 8 || m % 4 != 0)) continue;

      Matrix projected(m,m);
      for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = i; j < m; ++j)
          projected(i,j) = projected(j,i) = .5*(dot(q[i],aq[j])+dot(q[j],aq[i]));
      const auto eig = symmetric_eigen(projected);
      auto v = combination(q,eig.vectors,0);
      const double length = norm(v);
      for (auto& x : v) x /= length;
      // Reapply the original operator: do not rely solely on a projected
      // residual, which can hide accumulated orthogonality errors.
      const auto av = checked_apply(v);
      const double value = dot(v,av);
      Vector residual(size);
      for (std::size_t i = 0; i < size; ++i) residual[i] = av[i]-value*v[i];
      const double relative_residual = norm(residual)/std::max(1.0,std::abs(value));
      if (relative_residual <= options.tolerance)
        return {.value=value, .vector=std::move(v),
                .relative_residual=relative_residual, .products=products};

      if (m == capacity || exhausted) {
        Basis restarted;
        const auto retained = std::min<std::size_t>(4,std::max<std::size_t>(1,capacity/2));
        for (std::size_t axis = 0; axis < std::min(retained,m); ++axis)
          append_orthogonal(restarted,combination(q,eig.vectors,axis));
        append_orthogonal(restarted,std::move(residual));
        if (restarted.size() < capacity) append_orthogonal(restarted,random_vector());
        q = std::move(restarted);
        aq.clear();
        break;
      }
    }
  }
}
}
