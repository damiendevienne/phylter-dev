#pragma once
#include "phylter2/distatis.hpp"
#include "phylter2/preparation.hpp"
#include <functional>
#include <utility>
namespace phylter2 {
enum class WrNormalization { row, column, none };
struct AnalysisOptions {
  PreparationOptions preparation;
  double k{3}, k2{3}, stop{1e-5};
  bool islands{true}, initial_only{false};
  WrNormalization wr_normalization{WrNormalization::row};
  RvOptions rv{RvMethod::matrix_free, false};
};
using Cell = std::pair<std::size_t, std::size_t>; // gene, species
using StateObserver = std::function<void(std::size_t, const DistatisResult&, const Matrix&)>;
enum class OutlierMode { cells, whole_genes };
struct ProgressEvent {
  OutlierMode mode{OutlierMode::cells};
  std::size_t new_cells{};
  std::size_t new_genes{};
  double previous_quality{};
  double candidate_quality{};
  bool has_candidate{false};
  bool accepted{false};
};
using ProgressObserver = std::function<void(const ProgressEvent&)>;
struct AnalysisResult {
  std::vector<std::string> genes, taxa, discarded_genes;
  std::vector<std::pair<std::string,std::string>> outliers, discarded;
  std::vector<std::string> complete_genes, complete_species;
  std::vector<Cell> cells;
  std::vector<double> scores;
  // Compact final state: alpha, compromise, quality and eigenvalues only.
  // Use an observer to inspect full coordinates, centered matrices and RV.
  DistatisResult final_state;
  Matrix wr;
};
AnalysisResult analyze(std::vector<GeneMatrix> input, const AnalysisOptions& options = {},
                       const StateObserver& observer = {},
                       const ProgressObserver& progress = {});
}
