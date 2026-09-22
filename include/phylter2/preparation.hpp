#pragma once

#include "phylter2/newick.hpp"

#include <string>
#include <vector>

namespace phylter2 {

enum class Normalization { median, mean, none };

struct GeneMatrix {
  std::string gene;
  NamedDistanceMatrix distances;
};

struct PreparationOptions {
  Normalization normalization{Normalization::median};
  double normalization_cutoff{1e-3};
};

struct PreparationResult {
  std::vector<GeneMatrix> matrices;
  std::vector<std::string> discarded_genes;
  std::vector<std::string> taxa;
};
void validate_gene_input(const GeneMatrix& matrix);

// Imputes missing taxa to the ordered union of input taxa, then applies the
// exact R-package scale convention (median/mean over all matrix entries,
// including diagonal and mirrored distances). Never-co-occurring pairs use the
// R implementation's legacy nearest-neighbor fallback (see README).
[[nodiscard]] PreparationResult prepare_matrices(std::vector<GeneMatrix> matrices,
                                                 const PreparationOptions& options = {});

}  // namespace phylter2
