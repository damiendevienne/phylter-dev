#include "phylter2/preparation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace phylter2 {
namespace {

void validate_gene_matrix(const GeneMatrix& gene_matrix) {
  if (gene_matrix.gene.empty()) throw std::invalid_argument("gene identifiers must be non-empty");
  const NamedDistanceMatrix& distances = gene_matrix.distances;
  if (!distances.values.square() || distances.values.rows() != distances.taxa.size() ||
      distances.taxa.size() < 3) {
    throw std::invalid_argument("each gene needs a square distance matrix for at least three taxa");
  }
  std::unordered_set<std::string> labels;
  for (const std::string& taxon : distances.taxa) {
    if (taxon.empty() || !labels.insert(taxon).second) {
      throw std::invalid_argument("taxon labels must be unique and non-empty within each gene");
    }
  }
  for (std::size_t row = 0; row < distances.values.rows(); ++row) {
    if (std::abs(distances.values(row, row)) > 1e-12) {
      throw std::invalid_argument("distance-matrix diagonals must be zero");
    }
    for (std::size_t col = 0; col < distances.values.cols(); ++col) {
      const double value = distances.values(row, col);
      if (!std::isfinite(value) || value < 0.0 ||
          std::abs(value - distances.values(col, row)) > 1e-10) {
        throw std::invalid_argument("distance matrices must be finite, non-negative and symmetric");
      }
    }
  }
}

double median(std::vector<double> values) {
  if (values.empty()) throw std::invalid_argument("median of an empty matrix");
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle), values.end());
  const double upper = values[middle];
  if (values.size() % 2 != 0) return upper;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle - 1), values.end());
  return 0.5 * (upper + values[middle - 1]);
}

double scale(const Matrix& matrix, const Normalization normalization) {
  if (normalization == Normalization::none) return 1.0;
  if (normalization == Normalization::mean) {
    return std::accumulate(matrix.data().begin(), matrix.data().end(), 0.0) /
           static_cast<double>(matrix.data().size());
  }
  return median(std::vector<double>(matrix.data().begin(), matrix.data().end()));
}

}  // namespace

void validate_gene_input(const GeneMatrix& matrix) {
  validate_gene_matrix(matrix);
}

PreparationResult prepare_matrices(std::vector<GeneMatrix> matrices, const PreparationOptions& options) {
  if (matrices.size() < 2) throw std::invalid_argument("at least two genes are required");
  if (!std::isfinite(options.normalization_cutoff) || options.normalization_cutoff < 0.0) {
    throw std::invalid_argument("normalization cutoff must be a finite non-negative value");
  }

  std::unordered_set<std::string> gene_ids;
  std::vector<std::string> taxa;
  std::unordered_map<std::string, std::size_t> taxon_index;
  for (const GeneMatrix& matrix : matrices) {
    validate_gene_matrix(matrix);
    if (!gene_ids.insert(matrix.gene).second) throw std::invalid_argument("gene identifiers must be unique");
    for (const std::string& taxon : matrix.distances.taxa) {
      if (!taxon_index.contains(taxon)) {
        taxon_index.emplace(taxon, taxa.size());
        taxa.push_back(taxon);
      }
    }
  }

  const std::size_t species = taxa.size();
  Matrix total(species, species);
  Matrix count(species, species);
  for (const GeneMatrix& gene_matrix : matrices) {
    const NamedDistanceMatrix& input = gene_matrix.distances;
    for (std::size_t row = 0; row < input.taxa.size(); ++row) {
      const std::size_t expanded_row = taxon_index.at(input.taxa[row]);
      for (std::size_t col = 0; col < input.taxa.size(); ++col) {
        const std::size_t expanded_col = taxon_index.at(input.taxa[col]);
        total(expanded_row, expanded_col) += input.values(row, col);
        count(expanded_row, expanded_col) += 1.0;
      }
    }
  }

  Matrix average(species, species);
  for (std::size_t row = 0; row < species; ++row) {
    for (std::size_t col = 0; col < species; ++col) {
      average(row, col) = total(row, col) / count(row, col);
    }
  }

  // Compatibility with R impMean: both neighbor searches use the unmodified
  // mean matrix. Its mean(a,b) call treats b as 'trim', returning scalar a.
  // Consequently these fallback distances can be asymmetric. Do not silently
  // replace this with (a+b)/2 in an exact port.
  const Matrix& observed_average = average;
  std::vector<std::pair<std::size_t,double>> corrections;
  auto neighbor_distance = [&](std::size_t from, std::size_t toward) {
    std::vector<std::size_t> neighbors;
    for (std::size_t i = 0; i < species; ++i)
      if (std::isfinite(observed_average(toward,i))) neighbors.push_back(i);
    std::stable_sort(neighbors.begin(),neighbors.end(),[&](auto a, auto b) {
      return observed_average(toward,a) < observed_average(toward,b);
    });
    for (std::size_t i = 1; i < neighbors.size(); ++i)
      if (std::isfinite(observed_average(from,neighbors[i])))
        return observed_average(from,neighbors[i]);
    throw std::runtime_error("cannot impute a never-co-occurring taxon pair: no shared observed neighbor");
  };
  for (std::size_t row = 0; row < species; ++row)
    for (std::size_t col = 0; col < species; ++col)
      if (!std::isfinite(average(row,col))) {
        const double a = neighbor_distance(row,col);
        (void)neighbor_distance(col,row);
        corrections.emplace_back(row*species+col,a);
      }
  for (const auto& [index,value] : corrections) average.data()[index] = value;

  std::vector<GeneMatrix> prepared;
  std::vector<std::string> discarded;
  prepared.reserve(matrices.size());
  for (std::size_t gene = 0; gene < matrices.size(); ++gene) {
    // Expand one gene at a time: avoid an additional K full-sized matrices.
    Matrix matrix = average;
    auto& input = matrices[gene].distances;
    for (std::size_t row = 0; row < input.taxa.size(); ++row) {
      for (std::size_t col = 0; col < input.taxa.size(); ++col) {
        matrix(taxon_index.at(input.taxa[row]),taxon_index.at(input.taxa[col])) = input.values(row,col);
      }
    }
    input.values = Matrix();
    const double factor = scale(matrix, options.normalization);
    if (options.normalization != Normalization::none && factor <= options.normalization_cutoff) {
      discarded.push_back(matrices[gene].gene);
      continue;
    }
    if (!std::isfinite(factor) || factor == 0.0) {
      throw std::runtime_error("cannot normalize a matrix with zero or non-finite scale");
    }
    if (options.normalization != Normalization::none) {
      for (double& value : matrix.data()) value /= factor;
    }
    prepared.push_back({.gene = matrices[gene].gene,
                        .distances = {.values = std::move(matrix), .taxa = taxa}});
  }
  if (prepared.size() < 2) throw std::runtime_error("fewer than two genes remain after preparation");
  return {.matrices = std::move(prepared), .discarded_genes = std::move(discarded), .taxa = std::move(taxa)};
}

}  // namespace phylter2
