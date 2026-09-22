#include "phylter2/analysis.hpp"
#include "phylter2/statistics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
namespace phylter2 {
namespace {
std::vector<Cell> detect(const Matrix& wr, const DistatisResult& state,
                         const AnalysisOptions& options, bool whole_genes) {
  const auto n = wr.rows(), g = wr.cols();
  std::vector<Cell> result;
  if (whole_genes) {
    auto negative = state.alpha;
    for (auto& value : negative) value = -value;
    const double threshold = outlier_threshold(negative, options.k2);
    for (std::size_t gene = 0; gene < g; ++gene) if (negative[gene] > threshold)
      for (std::size_t species = 0; species < n; ++species) result.emplace_back(gene, species);
    return result;
  }
  const auto order = cluster_order(state.compromise);
  Matrix normalized = wr;
  if (options.wr_normalization != WrNormalization::none) {
    const bool row = options.wr_normalization == WrNormalization::row;
    for (std::size_t a = 0; a < (row ? n : g); ++a) {
      std::vector<double> values;
      for (std::size_t b = 0; b < (row ? g : n); ++b)
        values.push_back(row ? wr(a,b) : wr(b,a));
      const double median = quantile(std::move(values), .5);
      if (median == 0) throw std::runtime_error("WR normalization has a zero median; try --normalize-by none");
      for (std::size_t b = 0; b < (row ? g : n); ++b)
        (row ? normalized(a,b) : normalized(b,a)) /= median;
    }
  }
  std::vector<double> values;
  values.reserve(n*g);
  for (std::size_t gene = 0; gene < g; ++gene)
    for (auto species : order) values.push_back(normalized(species,gene));
  const double threshold = outlier_threshold(values, options.k);
  for (std::size_t gene = 0; gene < g; ++gene) {
    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < n;) {
      if (normalized(order[i],gene) <= threshold) { ++i; continue; }
      std::vector<std::size_t> group;
      do { group.push_back(order[i++]); }
      while (options.islands && i < n && normalized(order[i],gene) > threshold);
      groups.push_back(std::move(group));
    }
    if (options.islands)
      std::stable_partition(groups.begin(), groups.end(), [](const auto& group) { return group.size() > 1; });
    for (const auto& group : groups) {
      double maximum = -std::numeric_limits<double>::infinity();
      for (auto species : group) maximum = std::max(maximum, normalized(species,gene));
      for (auto species : group)
        if (normalized(species,gene) == maximum) result.emplace_back(gene,species);
    }
  }
  return result;
}
// Large diagnostic objects are transient in the CLI. Keep only what is needed
// for the next iteration; observers can serialize the full state first.
void compact(DistatisResult& state) {
  state.centered.clear();
  state.partial_coordinates.clear();
  state.coordinates = Matrix();
  state.rv = Matrix();
}
}
AnalysisResult analyze(std::vector<GeneMatrix> input, const AnalysisOptions& options,
                       const StateObserver& observer) {
  for (double value : {options.k, options.k2, options.stop})
    if (!std::isfinite(value) || value < 0) throw std::invalid_argument("k, k2 and stop must be finite and non-negative");
  std::unordered_map<std::string,std::vector<std::string>> original;
  for (const auto& gene : input) original.emplace(gene.gene, gene.distances.taxa);
  auto prepared = prepare_matrices(std::move(input), options.preparation);
  AnalysisResult result;
  result.taxa = prepared.taxa;
  result.discarded_genes = prepared.discarded_genes;
  for (const auto& gene : result.discarded_genes)
    for (const auto& species : original.at(gene)) result.discarded.emplace_back(gene,species);
  std::vector<Matrix> matrices;
  for (auto& gene : prepared.matrices) {
    result.genes.push_back(gene.gene);
    matrices.push_back(std::move(gene.distances.values));
  }
  auto state = distatis(matrices,0,-std::numeric_limits<double>::infinity(),options.rv);
  auto wr = dist2wr(state);
  result.scores.push_back(state.quality);
  if (observer) observer(0,state,wr);
  compact(state);
  const auto n = result.taxa.size(), g = matrices.size();
  std::vector<bool> known(n*g, false);
  bool whole_genes = false;
  while (!options.initial_only) {
    auto detected = detect(wr,state,options,whole_genes);
    std::vector<Cell> new_cells;
    for (auto cell : detected) if (!known[cell.first*n+cell.second]) new_cells.push_back(cell);
    if (new_cells.empty()) {
      if (whole_genes) break;
      whole_genes = true;
      continue;
    }
    auto cells = result.cells;
    cells.insert(cells.end(), new_cells.begin(), new_cells.end());
    auto proposal = matrices;
    // R replaces all selected rows first, then all columns, within each gene.
    for (auto [gene,species] : cells) for (std::size_t j = 0; j < n; ++j)
      proposal[gene](species,j) = state.compromise(species,j);
    for (auto [gene,species] : cells) for (std::size_t j = 0; j < n; ++j)
      proposal[gene](j,species) = state.compromise(j,species);
    // Preserve R's subtraction test exactly at the acceptance boundary.
    auto candidate = distatis(proposal, 0, state.quality + options.stop -
      4*std::numeric_limits<double>::epsilon(),options.rv);
    if (candidate.quality - state.quality < options.stop) {
      if (whole_genes) break;
      whole_genes = true;
      continue;
    }
    state = std::move(candidate);
    matrices = std::move(proposal);
    result.cells = std::move(cells);
    for (auto [gene,species] : new_cells) known[gene*n+species] = true;
    wr = dist2wr(state);
    result.scores.push_back(state.quality);
    if (observer) observer(result.scores.size()-1,state,wr);
    compact(state);
    whole_genes = false;
  }
  std::unordered_map<std::string,std::size_t> species_total, species_out;
  std::vector<std::unordered_set<std::string>> presence;
  for (const auto& gene : result.genes) {
    const auto& names = original.at(gene);
    presence.emplace_back(names.begin(),names.end());
    for (const auto& species : names) ++species_total[species];
  }
  std::vector<std::size_t> gene_out(g);
  for (auto [gene,species] : result.cells) if (presence[gene].contains(result.taxa[species])) {
    result.outliers.emplace_back(result.genes[gene],result.taxa[species]);
    ++gene_out[gene]; ++species_out[result.taxa[species]];
  }
  for (std::size_t gene = 0; gene < g; ++gene)
    if (gene_out[gene] == presence[gene].size()) result.complete_genes.push_back(result.genes[gene]);
  for (const auto& species : result.taxa)
    if (species_total[species] > 0 && species_out[species] == species_total[species])
      result.complete_species.push_back(species);
  std::sort(result.complete_genes.begin(),result.complete_genes.end());
  std::sort(result.complete_species.begin(),result.complete_species.end());
  result.final_state = std::move(state);
  result.wr = std::move(wr);
  return result;
}
}
