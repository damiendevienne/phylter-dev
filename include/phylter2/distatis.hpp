#pragma once

#include "phylter2/matrix.hpp"
#include "phylter2/rv.hpp"

#include <cstddef>
#include <limits>
#include <vector>

namespace phylter2 {

struct DistatisResult {
  Matrix coordinates;
  std::vector<Matrix> partial_coordinates;
  std::vector<double> alpha;
  std::vector<double> eigenvalues;
  Matrix rv;
  Matrix compromise;
  double quality{0.0};
  std::size_t rv_products{0};
  double rv_relative_residual{0};
  std::vector<Matrix> centered;
};

// `factors == 0` implements the current R package's automatic broken-stick
// choice. A positive value retains that many axes.
[[nodiscard]] DistatisResult distatis(const std::vector<Matrix>& matrices,
                                      std::size_t factors = 0,
                                      double minimum_quality = -std::numeric_limits<double>::infinity(),
                                      const RvOptions& rv_options = {});

[[nodiscard]] Matrix dist2wr(const DistatisResult& result);

}  // namespace phylter2
