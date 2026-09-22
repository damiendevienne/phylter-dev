#include "phylter2/analysis.hpp"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#ifdef _WIN32
#include <io.h>
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
namespace fs = std::filesystem;
using namespace phylter2;
namespace {
constexpr std::string_view version = "0.1.0-dev";
using Clock = std::chrono::steady_clock;
bool terminal(const int descriptor) {
#ifdef _WIN32
  return _isatty(descriptor) != 0;
#else
  return isatty(descriptor) != 0;
#endif
}
void enable_ansi(const unsigned long stream) {
#ifdef _WIN32
  const HANDLE handle = GetStdHandle(static_cast<DWORD>(stream));
  DWORD mode = 0;
  if (handle != INVALID_HANDLE_VALUE && GetConsoleMode(handle,&mode))
    SetConsoleMode(handle,mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#else
  (void)stream;
#endif
}
struct Style {
  bool color{terminal(
#ifdef _WIN32
    _fileno(stdout)
#else
    STDOUT_FILENO
#endif
  ) && std::getenv("NO_COLOR") == nullptr};
  std::string green(const std::string& text) const {
    return color ? "\033[32m"+text+"\033[0m" : text;
  }
  std::string yellow(const std::string& text) const {
    return color ? "\033[33m"+text+"\033[0m" : text;
  }
  std::string bold(const std::string& text) const {
    return color ? "\033[1m"+text+"\033[0m" : text;
  }
};
class Spinner {
 public:
  Spinner() : enabled_(terminal(
#ifdef _WIN32
    _fileno(stderr)
#else
    STDERR_FILENO
#endif
  )) {}
  ~Spinner() { stop(); }
  void start() {
    if (!enabled_ || running_) return;
    running_ = true;
    worker_ = std::thread([this] {
      constexpr std::string_view frames = "|/-\\";
      std::size_t frame = 0;
      while (running_) {
        std::cerr << "\r  " << frames[frame++%frames.size()] << " Computing..." << std::flush;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
    });
  }
  void stop() {
    if (!running_) return;
    running_ = false;
    if (worker_.joinable()) worker_.join();
    if (enabled_) std::cerr << "\r\033[2K" << std::flush;
  }
 private:
  bool enabled_{};
  std::atomic<bool> running_{false};
  std::thread worker_;
};
std::uint64_t available_memory() {
#ifdef _WIN32
  MEMORYSTATUSEX status{};
  status.dwLength = sizeof(status);
  return GlobalMemoryStatusEx(&status) ? status.ullAvailPhys : 0;
#elif defined(__APPLE__)
  vm_statistics64_data_t statistics{};
  mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
  vm_size_t page_size = 0;
  if (host_page_size(mach_host_self(),&page_size) != KERN_SUCCESS ||
      host_statistics64(mach_host_self(),HOST_VM_INFO64,
                        reinterpret_cast<host_info64_t>(&statistics),&count) != KERN_SUCCESS) return 0;
  return static_cast<std::uint64_t>(statistics.free_count+statistics.inactive_count)*page_size;
#else
  std::ifstream in("/proc/meminfo");
  std::string key,unit;
  std::uint64_t kib = 0;
  while (in >> key >> kib >> unit)
    if (key == "MemAvailable:") return kib*1024;
  return 0;
#endif
}
std::uint64_t peak_memory() {
#ifdef _WIN32
  PROCESS_MEMORY_COUNTERS counters{};
  return GetProcessMemoryInfo(GetCurrentProcess(),&counters,sizeof(counters))
    ? static_cast<std::uint64_t>(counters.PeakWorkingSetSize) : 0;
#else
  rusage usage{};
  if (getrusage(RUSAGE_SELF,&usage) != 0) return 0;
#ifdef __APPLE__
  return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
  return static_cast<std::uint64_t>(usage.ru_maxrss)*1024;
#endif
#endif
}
std::string memory_size(const long double bytes) {
  const char* units[] = {"B","KiB","MiB","GiB","TiB"};
  long double value = bytes;
  std::size_t unit = 0;
  while (value >= 1024 && unit < 4) { value /= 1024; ++unit; }
  std::ostringstream out;
  out << std::fixed << std::setprecision(value < 10 && unit ? 1 : 0)
      << static_cast<double>(value) << ' ' << units[unit];
  return out.str();
}
struct InputProfile {
  std::size_t trees{};
  std::size_t trees_with_lengths{};
  std::size_t trees_with_support{};
  std::size_t trees_with_usable_support{};
  std::vector<std::size_t> tips;
  std::vector<std::string> gene_names;
  std::unordered_set<std::string> taxa;
  std::unordered_map<std::string,std::size_t> tips_by_gene;
};
double median(std::vector<std::size_t> values) {
  std::sort(values.begin(),values.end());
  const auto middle = values.size()/2;
  if (values.size()%2) return static_cast<double>(values[middle]);
  return .5*static_cast<double>(values[middle-1]+values[middle]);
}
std::string decimal(const double value, const int precision = 6) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(precision) << value;
  return out.str();
}
std::string integer(const std::size_t value) {
  std::string result = std::to_string(value);
  for (std::ptrdiff_t position = static_cast<std::ptrdiff_t>(result.size())-3;
       position > 0; position -= 3)
    result.insert(static_cast<std::size_t>(position),",");
  return result;
}
void print_profile(const InputProfile& profile, const fs::path& input, const bool matrices,
                   const bool nodal, const double support_cutoff, const bool check,
                   const Style& style) {
  const auto [minimum,maximum] = std::minmax_element(profile.tips.begin(),profile.tips.end());
  const double coverage = 100.0*median(profile.tips)/static_cast<double>(profile.taxa.size());
  std::cout << style.green("Phylter "+std::string(version)) << "\n\n"
            << style.bold("Dataset") << '\n'
            << "  Input:             " << fs::absolute(input).string() << '\n'
            << "  Genes:             " << profile.tips.size() << '\n'
            << "  Taxa:              " << profile.taxa.size() << " unique\n"
            << "  Tips per gene:     median " << decimal(median(profile.tips),1)
            << ", range " << *minimum << '-' << *maximum
            << " (" << decimal(coverage,1) << "% median coverage)\n";
  if (matrices) {
    std::cout << "  Input type:        precomputed distance matrices\n";
  } else {
    std::cout << "  Trees:             " << profile.trees << '\n'
              << "  Branch lengths:    " << profile.trees_with_lengths << '/' << profile.trees
              << " trees\n"
              << "  Node support:      " << profile.trees_with_support << '/' << profile.trees
              << " trees";
    if (check) {
      std::cout << "\n\n" << style.bold("Available options") << '\n'
                << "  Nodal distance:    available (does not require branch lengths)\n"
                << "  Patristic distance: ";
      if (profile.trees_with_lengths == profile.trees) std::cout << "available\n";
      else std::cout << "unavailable; " << profile.trees-profile.trees_with_lengths
                     << " tree(s) have missing branch lengths\n";
      std::cout << "  Support collapse:  ";
      if (profile.trees_with_usable_support == profile.trees) std::cout << "available\n";
      else if (profile.trees_with_support == 0)
        std::cout << "unavailable; no internal-node support values found\n";
      else
        std::cout << "unavailable; numeric support is incomplete or invalid in "
                  << profile.trees-profile.trees_with_usable_support << " tree(s)\n";
      std::cout << '\n' << style.bold("Gene identifiers") << "\n  ";
      const std::size_t shown = std::min<std::size_t>(5,profile.gene_names.size());
      for (std::size_t i = 0; i < shown; ++i) {
        if (i) std::cout << ", ";
        std::cout << profile.gene_names[i];
      }
      if (shown < profile.gene_names.size())
        std::cout << ", ... (" << profile.gene_names.size() << " total)";
      std::cout << "\n  IDs use the filename stem; multi-tree files add :1, :2, and so on.\n";
    } else {
      if (support_cutoff == 0) std::cout << " (filtering disabled)\n";
      else std::cout << " (collapse below " << support_cutoff << ")\n";
      std::cout << "  Distance metric:   "
                << (nodal ? "nodal (counts edges; branch lengths are ignored)"
                          : "patristic (uses branch lengths)") << '\n';
    }
  }
  std::cout << '\n';
}
void help() {
  std::cout << R"(Phylter: standalone phylogenomic outlier filtering
Usage:
  phylter --trees FILE_OR_DIRECTORY --out PREFIX [options]
  phylter --matrices DIRECTORY --out PREFIX [options]
  phylter --check --trees FILE_OR_DIRECTORY

Input: Newick file with one or more trees, or directory of Newick files.
Matrices: tab-separated *.tsv files, first row/column are taxon names.
Gene names: file stem; stem:1, stem:2, ... for multiple trees in one file.

Options:
  --distance patristic|nodal  Tree distances (default: patristic)
  --norm median|mean|none     Distance normalization (default: median)
  --norm-cutoff NUMBER       Discard uninformative genes (default: 0.001)
  --support-cutoff NUMBER    Collapse internal branches below this support (default: 0)
  --k NUMBER                 Cell outlier threshold (default: 3)
  --k2 NUMBER                Whole gene threshold (default: same as k)
  --stop NUMBER              Minimum quality gain (default: 0.00001)
  --normalize-by row|col|none WR normalization (default: row)
  --no-islands               Remove every flagged cell
  --initial-only             Compute the initial state only
  --diagnostics DIRECTORY    Write each accepted state's matrices (large)
  --check                    Inspect tree structure and report available analyses
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
  const auto program_started = Clock::now();
#ifdef _WIN32
  enable_ansi(STD_OUTPUT_HANDLE);
  enable_ansi(STD_ERROR_HANDLE);
#endif
  const Style style;
  try {
    if (argc == 1) { help(); return 0; }
    for (int i = 1; i < argc; ++i)
      if (std::string_view(argv[i]) == "--help" || std::string_view(argv[i]) == "-h") { help(); return 0; }
    const bool check_requested = std::any_of(argv+1,argv+argc,[](const char* argument) {
      return std::string_view(argument) == "--check";
    });
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
      std::cout << "phylter " << version << '\n';
      return 0;
    }
    AnalysisOptions options;
    fs::path input, prefix, diagnostics;
    bool matrices = false, nodal = false, force = false, check = check_requested, k2_set = false;
    double support_cutoff = 0;
    int threads = 1;
    std::unordered_set<std::string> seen;
    for (int i = 1; i < argc; ++i) {
      const std::string key = argv[i];
      if (!seen.insert(key).second) throw std::invalid_argument("repeated option: " + key);
      if (check && key != "--check" && key != "--trees")
        throw std::invalid_argument("--check only accepts --trees");
      if (key == "--no-islands") { options.islands = false; continue; }
      if (key == "--initial-only") { options.initial_only = true; continue; }
      if (key == "--check") continue;
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
    if (!check && prefix.empty()) throw std::invalid_argument("--out is required unless --check is used");
    if (matrices && nodal) throw std::invalid_argument("--distance applies only to trees");
    std::vector<GeneMatrix> genes;
    InputProfile profile;
    for (const auto& path : input_files(input,matrices)) {
      try {
        if (matrices) genes.push_back({path.stem().string(),read_matrix(path)});
        else {
          auto trees = parse_newick_set(read_file(path));
          for (std::size_t t = 0; t < trees.size(); ++t) {
            ++profile.trees;
            const auto& tree = trees[t];
            bool complete_lengths = true;
            bool has_support = false;
            bool usable_support = true;
            for (std::size_t node_index = 0; node_index < tree.nodes.size(); ++node_index) {
              const auto& node = tree.nodes[node_index];
              if (node_index != tree.root) {
                if (!node.has_length) complete_lengths = false;
                else if (!std::isfinite(node.length_to_parent) || node.length_to_parent < 0)
                  throw std::invalid_argument("branch lengths must be finite and non-negative");
              }
              if (!node.children.empty()) {
                if (!node.label.empty()) has_support = true;
                std::size_t end = 0;
                try {
                  const double value = std::stod(node.label,&end);
                  if (end != node.label.size() || !std::isfinite(value)) usable_support = false;
                } catch (const std::exception&) {
                  usable_support = false;
                }
              }
            }
            if (complete_lengths) ++profile.trees_with_lengths;
            if (has_support) ++profile.trees_with_support;
            if (has_support && usable_support) ++profile.trees_with_usable_support;
            const auto name = path.stem().string() + (trees.size() > 1 ? ":" + std::to_string(t+1) : "");
            genes.push_back({name,patristic_distances(trees[t],check ? true : nodal,
                                                      check ? 0.0 : support_cutoff)});
          }
        }
      } catch (const std::exception& error) {
        throw std::runtime_error(path.string() + ": " + error.what());
      }
    }
    std::unordered_set<std::string> ids;
    for (const auto& gene : genes) {
      label_check(gene.gene);
      try {
        validate_gene_input(gene);
      } catch (const std::exception& error) {
        throw std::runtime_error("gene '" + gene.gene + "': " + error.what());
      }
      if (!ids.insert(gene.gene).second) throw std::invalid_argument("duplicate gene identifier");
      profile.tips.push_back(gene.distances.taxa.size());
      profile.gene_names.push_back(gene.gene);
      profile.tips_by_gene.emplace(gene.gene,gene.distances.taxa.size());
      for (const auto& taxon : gene.distances.taxa) {
        label_check(taxon);
        profile.taxa.insert(taxon);
      }
    }
    print_profile(profile,input,matrices,nodal,support_cutoff,check,style);
    if (check) {
      std::cout << style.green("The tree dataset is well formed.") << '\n';
      return 0;
    }
    for (const std::string suffix : {".outliers.tsv",".discarded.tsv",".scores.tsv",".summary.txt"})
      if (!force && fs::exists(prefix.string()+suffix)) throw std::runtime_error("output exists (use --force): " + prefix.string()+suffix);
    if (!prefix.parent_path().empty()) fs::create_directories(prefix.parent_path());
    if (!diagnostics.empty()) {
      if (fs::exists(diagnostics) && (!fs::is_directory(diagnostics) || !fs::is_empty(diagnostics)))
        throw std::runtime_error("diagnostics directory must be empty");
      fs::create_directories(diagnostics);
    }
    const char* normalization = options.preparation.normalization == Normalization::median ? "median" :
      options.preparation.normalization == Normalization::mean ? "mean" : "none";
    const long double gene_count = static_cast<long double>(profile.tips.size());
    const long double taxon_count = static_cast<long double>(profile.taxa.size());
    long double estimated_memory = 6.0L*gene_count*taxon_count*taxon_count*sizeof(double);
    estimated_memory += 3.0L*taxon_count*taxon_count*sizeof(double);
    if (options.rv.method == RvMethod::dense || !diagnostics.empty())
      estimated_memory += 2.0L*gene_count*gene_count*sizeof(double);
    else
      estimated_memory += 96.0L*gene_count*sizeof(double);
    const std::uint64_t memory_available = available_memory();
    const long double memory_ratio = memory_available == 0 ? 0 :
      estimated_memory/static_cast<long double>(memory_available);
    std::cout << style.bold("Settings (configurable; run phylter --help for all options)") << '\n'
              << "  Detection:         k=" << options.k << ", k2=" << options.k2
              << ", " << (options.islands ? "island maxima" : "all flagged cells") << '\n'
              << "  Distance scaling:  ";
    if (options.preparation.normalization == Normalization::none)
      std::cout << "disabled (--norm none)\n"
                << "  Gene filtering:    disabled when distance scaling is disabled\n";
    else
      std::cout << "divide each gene matrix by its " << normalization
                << " distance (--norm " << normalization << ")\n"
                << "  Gene filtering:    discard genes with " << normalization << " distance <= "
                << options.preparation.normalization_cutoff << " (--norm-cutoff)\n";
    std::cout << "  Acceptance rule:   apply removals when quality gain is at least "
              << decimal(options.stop) << " (--stop)\n"
              << "  RV solver:         "
              << (options.rv.method == RvMethod::matrix_free ? "matrix-free" : "dense")
              << ", " << threads << " thread(s)\n"
              << "  Estimated memory:  ~" << memory_size(estimated_memory)
              << " (rough upper estimate)";
    if (memory_available)
      std::cout << ", " << decimal(100.0*static_cast<double>(memory_ratio),1)
                << "% of currently available RAM";
    std::cout << '\n';
    if (memory_ratio >= .80L)
      std::cout << style.yellow("  WARNING: this run may exhaust available memory; use a machine with more RAM.") << '\n';
    else if (memory_ratio >= .50L)
      std::cout << style.yellow("  CAUTION: this run may put substantial pressure on available memory.") << '\n';
    std::cout << '\n';
    options.rv.materialize = !diagnostics.empty();
    Spinner spinner;
    auto observer = [&](std::size_t iteration, const DistatisResult& state, const Matrix& wr) {
      if (iteration == 0) {
        spinner.stop();
        std::cout << style.bold("Optimization") << "\n  Initial quality:   "
                  << decimal(state.quality) << '\n' << std::flush;
        spinner.start();
      }
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
    auto progress = [&](const ProgressEvent& event) {
      spinner.stop();
      const bool genes_mode = event.mode == OutlierMode::whole_genes;
      if (!event.has_candidate) {
        std::cout << "  No new " << (genes_mode ? "complete-gene" : "cell")
                  << " outliers | "
                  << (genes_mode ? "optimization complete" : "checking complete-gene outliers")
                  << ".\n" << std::flush;
        spinner.start();
        return;
      }
      std::cout << "  ";
      if (genes_mode)
        std::cout << integer(event.new_genes) << " candidate complete-gene outlier"
                  << (event.new_genes == 1 ? "" : "s") << " ("
                  << integer(event.new_cells) << " cells)";
      else
        std::cout << integer(event.new_cells) << " candidate cell outlier"
                  << (event.new_cells == 1 ? "" : "s");
      const double gain = event.candidate_quality-event.previous_quality;
      std::cout << " | quality gain " << (gain >= 0 ? "+" : "") << decimal(gain);
      if (event.accepted)
        std::cout << " | " << style.green("removal applied")
                  << " | new quality " << decimal(event.candidate_quality) << '\n';
      else
        std::cout << " (too small) | " << style.yellow("no removal")
                  << " | quality remains " << decimal(event.previous_quality) << '\n';
      if (!event.accepted)
        std::cout << "    " << (genes_mode ? "Optimization complete"
                                            : "Checking complete-gene outliers") << ".\n";
      std::cout << std::flush;
      spinner.start();
    };
    spinner.start();
    const auto result = analyze(std::move(genes), options, observer, progress);
    spinner.stop();
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
    const std::size_t analyzed_cells = std::accumulate(result.genes.begin(),result.genes.end(),std::size_t{},
      [&](const std::size_t total,const std::string& gene) { return total+profile.tips_by_gene.at(gene); });
    const double filtered = analyzed_cells == 0 ? 0.0 :
      100.0*static_cast<double>(result.outliers.size())/static_cast<double>(analyzed_cells);
    const double elapsed = std::chrono::duration<double>(Clock::now()-program_started).count();
    const std::uint64_t memory_peak = peak_memory();
    std::cout << '\n' << style.bold("Results") << '\n'
              << "  Genes analyzed:    " << result.genes.size();
    if (!result.discarded_genes.empty())
      std::cout << " (" << result.discarded_genes.size() << " discarded during preparation)";
    std::cout << "\n  Accepted cycles:   " << result.scores.size()-1
              << "\n  Final quality:     " << decimal(result.scores.back())
              << "\n  Quality gain:      " << decimal(100.0*(result.scores.back()-result.scores.front()),2)
              << '%'
              << "\n  Data loss:         " << decimal(filtered,2) << "% ("
              << result.outliers.size() << '/' << analyzed_cells << " outlier cells removed)"
              << "\n  Complete genes:    " << result.complete_genes.size()
              << "\n  Complete species:  " << result.complete_species.size()
              << "\n  Elapsed time:      " << decimal(elapsed,2) << " s";
    if (memory_peak) std::cout << "\n  Peak memory:       " << memory_size(memory_peak);
    std::cout << "\n\n" << style.bold("Output files") << '\n'
              << "  Summary:           " << fs::absolute(prefix.string()+".summary.txt").string() << '\n'
              << "  Outlier pairs:     " << fs::absolute(prefix.string()+".outliers.tsv").string() << '\n'
              << "  Discarded pairs:   " << fs::absolute(prefix.string()+".discarded.tsv").string() << '\n'
              << "  Quality scores:    " << fs::absolute(prefix.string()+".scores.tsv").string() << "\n\n"
              << "Note: rerun with --diagnostics DIR to write per-state matrices to a directory "
                 "(large output).\n";
    return 0;
  } catch (const std::bad_alloc&) {
    std::cerr << "phylter: not enough memory to complete the analysis; use a machine with more RAM "
                 "or reduce the number of genes/taxa\n";
    return 1;
  } catch (const std::exception& error) {
    std::cerr << "phylter: " << error.what() << '\n';
    return 1;
  }
}
