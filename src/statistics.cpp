#include "phylter2/statistics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
namespace phylter2 {
double quantile(std::vector<double> values, double probability) {
  if (values.empty() || !std::isfinite(probability) || probability < 0 || probability > 1)
    throw std::invalid_argument("invalid quantile input");
  for (double v : values) if (!std::isfinite(v))
    throw std::invalid_argument("quantile requires finite values");
  const double position = static_cast<double>(values.size() - 1) * probability;
  const auto lower = static_cast<std::size_t>(position);
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(lower), values.end());
  const double a = values[lower], fraction = position - static_cast<double>(lower);
  if (fraction == 0) return a;
  const double b = *std::min_element(values.begin() + static_cast<std::ptrdiff_t>(lower + 1), values.end());
  return (1 - fraction) * a + fraction * b; // R quantile type 7
}
double outlier_threshold(const std::vector<double>& values, double k) {
  const double q3 = quantile(values, .75);
  return q3 + k * std::exp(3 * medcouple(values)) * (q3 - quantile(values, .25)) + 1e-10;
}
// Complete linkage, with R hclust's right-neighbor and merge-order tie rules.
// stats/src/hclust.f (F. Murtagh and R contributors), GPL-2-or-later.
std::vector<std::size_t> cluster_order(const Matrix& distances) {
  const auto n = distances.rows();
  if (!distances.square() || n < 2) throw std::invalid_argument("clustering needs a square matrix");
  Matrix d = distances;
  for (std::size_t i = 0; i < n; ++i) for (std::size_t j = i + 1; j < n; ++j) {
    if (!std::isfinite(d(j,i))) throw std::invalid_argument("non-finite clustering distance");
    d(i,j) = d(j,i); // R as.dist uses the lower triangle
  }
  std::vector<bool> active(n, true);
  std::vector<std::size_t> neighbor(n), node(n);
  std::iota(node.begin(), node.end(), 0);
  std::vector<double> nearest(n, std::numeric_limits<double>::infinity());
  std::vector<std::pair<std::size_t,std::size_t>> merges;
  auto update = [&](std::size_t i) {
    nearest[i] = std::numeric_limits<double>::infinity();
    for (std::size_t j = i+1; j < n; ++j)
      if (active[j] && d(i,j) < nearest[i]) { nearest[i] = d(i,j); neighbor[i] = j; }
  };
  for (std::size_t i = 0; i < n; ++i) update(i);
  for (std::size_t step = 0; step < n-1; ++step) {
    std::size_t i = n;
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < n; ++k)
      if (active[k] && nearest[k] < best) { best = nearest[k]; i = k; }
    if (i == n) throw std::runtime_error("clustering could not find a finite merge");
    const auto j = neighbor[i];
    auto a = node[i], b = node[j];
    if ((a >= n && b < n) || (a >= n && b >= n && a > b)) std::swap(a,b);
    merges.emplace_back(a,b);
    node[i] = n + step;
    active[j] = false;
    for (std::size_t k = 0; k < n; ++k) if (active[k] && k != i)
      d(i,k) = d(k,i) = std::max(d(i,k), d(j,k));
    update(i);
    for (std::size_t k = 0; k < n; ++k)
      if (active[k] && (neighbor[k] == i || neighbor[k] == j)) update(k);
  }
  std::vector<std::size_t> result, stack{2*n-2};
  while (!stack.empty()) {
    const auto id = stack.back(); stack.pop_back();
    if (id < n) result.push_back(id);
    else { stack.push_back(merges[id-n].second); stack.push_back(merges[id-n].first); }
  }
  return result;
}
}
