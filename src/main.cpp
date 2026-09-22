#include "phylter2/analysis.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace fs = std::filesystem;
using namespace phylter2;
namespace {
constexpr std::string_view version = "0.1.0-dev";
void help() {
  std::cout << R"(PhylteR: standalone phylogenomic outlier filtering
Usage:
  phylter run --trees FILE_OR_DIRECTORY --out PREFIX [options]
  phylter run --matrices DIRECTORY --out PREFIX [options]
  phylter check --trees FILE_OR_DIRECTORY [--distance nodal]

Input: Newick file with one or more trees, or directory of Newick files.
Matrices: tab-separated *.tsv files, first row/column are taxon names.
Gene names: file stem; stem:1, stem:2, ... for multiple trees in one file.

Options:
  --distance patristic|nodal  Tree distances (default: patristic)
  --norm median|mean|none     Distance normalization (default: median)
  --norm-cutoff NUMBER       Discard uninformative genes (default: 0.001)
  --support-cutoff NUMBER    R-compatible node support filtering (default: 0)
  --k NUMBER                 Cell outlier threshold (default: 3)
  --k2 NUMBER                Whole gene threshold (default: same as k)
  --stop NUMBER              Minimum quality gain (default: 0.00001)
  --normalize-by row|col|none WR normalization (default: row)
  --no-islands               Remove every flagged cell
  --initial-only             Compute the initial state only
  --diagnostics DIRECTORY    Write each accepted state's matrices (large)
  --force                    Overwrite existing output files
  --threads INTEGER          RV worker threads (default: 1)
  --rv-method matrix-free|dense  RV solver (default: matrix-free)
  --help, --version

Outputs: PREFIX.outliers.tsv, .discarded.tsv, .scores.tsv, .summary.txt.
Compact output by default. This development version uses dense matrices.
)";
}
std::string read_file(const fs::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path.string());
  std::ostringstream out; out << in.rdbuf();
  if (in.bad()) throw std::runtime_error("error reading " + path.string());
  return out.str();
}
std::vector<fs::path> input_files(const fs::path& input, bool matrices) {
  if (fs::is_regular_file(input)) return {input};
  if (!fs::is_directory(input)) throw std::runtime_error("input does not exist: " + input.string());
  std::vector<fs::path> files;
  for (const auto& entry : fs::directory_iterator(input)) {
    if (!entry.is_regular_file()) continue;
    auto ext = entry.path().extension().string();
    std::transform(ext.begin(),ext.end(),ext.begin(),[](unsigned char c){return static_cast<char>(std::tolower(c));});
    if (matrices ? ext == ".tsv" : (ext == ".nwk" || ext == ".newick" ||
        ext == ".tre" || ext == ".tree" || ext == ".treefile")) files.push_back(entry.path());
  }
  std::sort(files.begin(),files.end());
  if (files.empty()) throw std::runtime_error("no supported input files");
  return files;
}
double number(const std::string& text) {
  std::size_t end;
  const double value = std::stod(text,&end);
  if (end != text.size()) throw std::invalid_argument("invalid number: " + text);
  return value;
}
std::vector<std::string> fields(const std::string& line) {
  std::vector<std::string> values; std::istringstream in(line); std::string value;
  while (std::getline(in,value,'\t')) {
    if (!value.empty() && value.back() == '\r') value.pop_back();
    values.push_back(value);
  }
  return values;
}
NamedDistanceMatrix read_matrix(const fs::path& path) {
  std::istringstream in(read_file(path)); std::string line;
  if (!std::getline(in,line)) throw std::runtime_error("empty matrix");
  auto names = fields(line);
  if (names.size() < 4) throw std::runtime_error("matrix header requires row label and at least three taxa");
  names.erase(names.begin());
  Matrix values(names.size(),names.size());
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (!std::getline(in,line)) throw std::runtime_error("missing matrix row");
    const auto row = fields(line);
    if (row.size() != names.size()+1 || row[0] != names[i])
      throw std::runtime_error("matrix row labels/order must match the header");
    for (std::size_t j = 0; j < names.size(); ++j) values(i,j) = number(row[j+1]);
  }
  while (std::getline(in,line)) if (!line.empty()) throw std::runtime_error("extra matrix rows");
  return {.values = std::move(values), .taxa = std::move(names)};
}
void label_check(const std::string& value) {
  if (value.empty() || value.find_first_of("\t\r\n") != std::string::npos)
    throw std::invalid_argument("identifiers cannot be empty or contain tabs/newlines");
}
std::ofstream output(const fs::path& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path.string());
  out.exceptions(std::ios::badbit | std::ios::failbit);
  out << std::setprecision(17);
  return out;
}
void write_matrix(const fs::path& path, const Matrix& matrix) {
  auto out = output(path);
  for (std::size_t i = 0; i < matrix.rows(); ++i) {
    for (std::size_t j = 0; j < matrix.cols(); ++j) {
      if (j) out << '\t';
      out << matrix(i,j);
    }
    out << '\n';
  }
}
void write_pairs(const fs::path& path, const std::vector<std::pair<std::string,std::string>>& pairs) {
  auto out = output(path); out << "gene\tspecies\n";
  for (const auto& [gene,species] : pairs) out << gene << '\t' << species << '\n';
}
}
int main(int argc, char* argv[]) {
  try {
    if (argc == 1) { help(); return 0; }
    for (int i = 1; i < argc; ++i)
      if (std::string_view(argv[i]) == "--help" || std::string_view(argv[i]) == "-h") { help(); return 0; }
    const std::string command = argv[1];
    if (command == "--version") { std::cout << "phylter " << version << '\n'; return 0; }
    if (command != "run" && command != "check") throw std::invalid_argument("unknown command: " + command);
    AnalysisOptions options;
    fs::path input, prefix, diagnostics;
    bool matrices = false, nodal = false, force = false, k2_set = false;
    double support_cutoff = 0;
    int threads = 1;
    std::unordered_set<std::string> seen;
    for (int i = 2; i < argc; ++i) {
      const std::string key = argv[i];
      if (!seen.insert(key).second) throw std::invalid_argument("repeated option: " + key);
      if (key == "--no-islands") { options.islands = false; continue; }
      if (key == "--initial-only") { options.initial_only = true; continue; }
      if (key == "--force") { force = true; continue; }
      if (i+1 == argc) throw std::invalid_argument("missing value for " + key);
      const std::string value = argv[++i];
      if (key == "--trees" || key == "--matrices") {
        if (!input.empty()) throw std::invalid_argument("choose either --trees or --matrices");
        input = value; matrices = key == "--matrices";
      } else if (key == "--out") prefix = value;
      else if (key == "--diagnostics") diagnostics = value;
      else if (key == "--k") options.k = number(value);
      else if (key == "--k2") { options.k2 = number(value); k2_set = true; }
      else if (key == "--stop") options.stop = number(value);
      else if (key == "--norm-cutoff") options.preparation.normalization_cutoff = number(value);
      else if (key == "--support-cutoff") support_cutoff = number(value);
      else if (key == "--rv-method") {
        if (value == "matrix-free") options.rv.method = RvMethod::matrix_free;
        else if (value == "dense") options.rv.method = RvMethod::dense;
        else throw std::invalid_argument("invalid RV method: " + value);
      }
      else if (key == "--threads") {
        std::size_t end;
        threads = std::stoi(value,&end);
        if (end != value.size() || threads < 1) throw std::invalid_argument("threads must be a positive integer");
      }
      else if (key == "--distance") {
        if (value != "patristic" && value != "nodal") throw std::invalid_argument("invalid distance");
        nodal = value == "nodal";
      } else if (key == "--norm") {
        if (value == "median") options.preparation.normalization = Normalization::median;
        else if (value == "mean") options.preparation.normalization = Normalization::mean;
        else if (value == "none") options.preparation.normalization = Normalization::none;
        else throw std::invalid_argument("invalid normalization");
      } else if (key == "--normalize-by") {
        if (value == "row") options.wr_normalization = WrNormalization::row;
        else if (value == "col") options.wr_normalization = WrNormalization::column;
        else if (value == "none") options.wr_normalization = WrNormalization::none;
        else throw std::invalid_argument("invalid WR normalization");
      } else throw std::invalid_argument("unknown option: " + key);
    }
    if (!k2_set) options.k2 = options.k;
    #ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
    #else
    if (threads != 1) throw std::runtime_error("this build has no OpenMP support; use --threads 1");
    #endif
    if (input.empty()) throw std::invalid_argument("--trees or --matrices is required");
    if (command == "run" && prefix.empty()) throw std::invalid_argument("--out is required");
    if (matrices && nodal) throw std::invalid_argument("--distance applies only to trees");
    std::vector<GeneMatrix> genes;
    for (const auto& path : input_files(input,matrices)) {
      if (matrices) genes.push_back({path.stem().string(),read_matrix(path)});
      else {
        auto trees = parse_newick_set(read_file(path));
        for (std::size_t t = 0; t < trees.size(); ++t) {
          const auto name = path.stem().string() + (trees.size() > 1 ? ":" + std::to_string(t+1) : "");
          genes.push_back({name,patristic_distances(trees[t],nodal,support_cutoff)});
        }
      }
    }
    std::unordered_set<std::string> ids;
    for (const auto& gene : genes) {
      label_check(gene.gene);
      validate_gene_input(gene);
      if (!ids.insert(gene.gene).second) throw std::invalid_argument("duplicate gene identifier");
      for (const auto& taxon : gene.distances.taxa) label_check(taxon);
    }
    if (command == "check") { std::cout << "Valid input: " << genes.size() << " genes\n"; return 0; }
    for (const std::string suffix : {".outliers.tsv",".discarded.tsv",".scores.tsv",".summary.txt"})
      if (!force && fs::exists(prefix.string()+suffix)) throw std::runtime_error("output exists (use --force): " + prefix.string()+suffix);
    if (!prefix.parent_path().empty()) fs::create_directories(prefix.parent_path());
    if (!diagnostics.empty()) {
      if (fs::exists(diagnostics) && (!fs::is_directory(diagnostics) || !fs::is_empty(diagnostics)))
        throw std::runtime_error("diagnostics directory must be empty");
      fs::create_directories(diagnostics);
    }
    options.rv.materialize = !diagnostics.empty();
    auto observer = [&](std::size_t iteration, const DistatisResult& state, const Matrix& wr) {
      std::cerr << "State " << iteration << ": quality " << std::setprecision(12) << state.quality << '\n';
      if (diagnostics.empty()) return;
      const auto base = diagnostics / ("state-" + std::to_string(iteration));
      write_matrix(base.string()+".wr.tsv",wr);
      write_matrix(base.string()+".rv.tsv",state.rv);
      {
        auto out = output(base.string()+".rv-solver.tsv");
        out << "products\trelative_residual\n" << state.rv_products << '\t';
        if (state.rv_products) out << state.rv_relative_residual;
        else out << "NA";
        out << '\n';
      }
      write_matrix(base.string()+".compromise.tsv",state.compromise);
      Matrix gram(state.coordinates.rows(),state.coordinates.rows());
      for (std::size_t i = 0; i < gram.rows(); ++i) for (std::size_t j = 0; j < gram.cols(); ++j)
        for (std::size_t a = 0; a < state.coordinates.cols(); ++a)
          gram(i,j) += state.coordinates(i,a)*state.coordinates(j,a);
      write_matrix(base.string()+".fgram.tsv",gram);
      auto out = output(base.string()+".weights.tsv");
      for (double a : state.alpha) out << a << '\n';
    };
    const auto result = analyze(std::move(genes), options, observer);
    write_pairs(prefix.string()+".outliers.tsv",result.outliers);
    write_pairs(prefix.string()+".discarded.tsv",result.discarded);
    {
      auto out = output(prefix.string()+".scores.tsv"); out << "state\tquality\n";
      for (std::size_t i = 0; i < result.scores.size(); ++i) out << i << '\t' << result.scores[i] << '\n';
    }
    {
      auto out = output(prefix.string()+".summary.txt");
      out << "phylter " << version << "\nGenes analyzed: " << result.genes.size()
          << "\nTaxa: " << result.taxa.size() << "\nDiscarded genes: " << result.discarded_genes.size()
          << "\nOutlier gene/species pairs: " << result.outliers.size()
          << "\nAccepted iterations: " << result.scores.size()-1
          << "\nInitial quality: " << result.scores.front()
          << "\nFinal quality: " << result.scores.back() << "\nComplete gene outliers:";
      for (const auto& gene : result.complete_genes) out << ' ' << gene;
      out << "\nComplete species outliers:";
      for (const auto& species : result.complete_species) out << ' ' << species;
      const char* norm = options.preparation.normalization == Normalization::median ? "median" :
        options.preparation.normalization == Normalization::mean ? "mean" : "none";
      const char* wrnorm = options.wr_normalization == WrNormalization::row ? "row" :
        options.wr_normalization == WrNormalization::column ? "col" : "none";
      out << "\nInput: " << fs::absolute(input).string()
          << "\nDistance: " << (matrices ? "provided matrices" : nodal ? "nodal" : "patristic")
          << "\nNormalization: " << norm << "\nNormalization cutoff: " << options.preparation.normalization_cutoff
          << "\nWR normalization: " << wrnorm << "\nk: " << options.k << "\nk2: " << options.k2
          << "\nMinimum gain: " << options.stop << "\nSupport cutoff: " << support_cutoff
          << "\nIslands: " << options.islands << "\nInitial only: " << options.initial_only
          << "\nRV threads: " << threads
          << "\nRV method: " << (options.rv.method == RvMethod::matrix_free ? "matrix-free" : "dense") << '\n';
    }
    if (!diagnostics.empty()) {
      auto cells = output(diagnostics/"cells.tsv");
      for (auto [gene,species] : result.cells) cells << gene+1 << '\t' << species+1 << '\n';
      auto taxa = output(diagnostics/"taxa.txt");
      for (const auto& taxon : result.taxa) taxa << taxon << '\n';
      auto names = output(diagnostics/"genes.txt");
      for (const auto& gene : result.genes) names << gene << '\n';
    }
    std::cout << result.outliers.size() << " outlier pairs; results: " << prefix.string() << ".*\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "phylter: " << error.what() << '\n';
    return 1;
  }
}
