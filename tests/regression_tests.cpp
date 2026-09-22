#include "phylter2/analysis.hpp"
#include "phylter2/statistics.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace phylter2;
namespace fs = std::filesystem;
namespace {
void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
std::ifstream input(const fs::path& path) {
  std::ifstream in(path);
  require(bool(in),"cannot read fixture " + path.string());
  return in;
}
std::vector<double> numbers(const fs::path& path) {
  auto in = input(path); std::vector<double> result; double value;
  while (in >> value) result.push_back(value);
  require(in.eof(),"invalid numeric fixture");
  return result;
}
void compare(const Matrix& actual, const fs::path& path) {
  const auto expected = numbers(path);
  require(expected.size() == actual.data().size(),"matrix size mismatch: " + path.string());
  for (std::size_t i = 0; i < expected.size(); ++i)
    require(std::isfinite(actual.data()[i]) && std::abs(actual.data()[i]-expected[i]) < 1e-8,
            "matrix value mismatch: " + path.string());
}
}
int main(int argc, char* argv[]) {
  try {
    require(argc==2,"fixture directory required");
    #ifdef _OPENMP
    omp_set_num_threads(1);
    #endif
    const fs::path root = argv[1];
    auto in = input(root/"carnivora.nwk"); std::ostringstream buffer; buffer << in.rdbuf();
    auto trees = parse_newick_set(buffer.str());
    std::vector<GeneMatrix> genes;
    for (std::size_t i = 0; i < trees.size(); ++i)
      genes.push_back({"carnivora:" + std::to_string(i+1),patristic_distances(trees[i])});
    const auto scores = numbers(root/"carnivora.scores.tsv");
    auto observer = [&](std::size_t iteration,const DistatisResult& state,const Matrix& wr) {
      require(iteration < scores.size(),"too many accepted states");
      require(std::abs(state.quality-scores[iteration]) < 1e-10,"quality trajectory mismatch");
      if (iteration!=0 && iteration!=scores.size()-1) return;
      const auto base = root / (iteration==0 ? "carnivora.initial" : "carnivora.final");
      compare(wr,base.string()+".wr.tsv");
      compare(state.rv,base.string()+".rv.tsv");
      compare(state.compromise,base.string()+".compromise.tsv");
      Matrix weights(state.alpha.size(),1);
      std::copy(state.alpha.begin(),state.alpha.end(),weights.data().begin());
      compare(weights,base.string()+".weights.tsv");
      Matrix gram(state.coordinates.rows(),state.coordinates.rows());
      for (std::size_t i = 0; i < gram.rows(); ++i) for (std::size_t j = 0; j < gram.cols(); ++j)
        for (std::size_t a = 0; a < state.coordinates.cols(); ++a)
          gram(i,j) += state.coordinates(i,a)*state.coordinates(j,a);
      compare(gram,base.string()+".fgram.tsv");
    };
    AnalysisOptions options;
    options.rv.materialize = true;
    const auto result = analyze(std::move(genes),options,observer);
    require(result.scores.size()==scores.size(),"wrong number of accepted states");
    auto pairs = input(root/"carnivora.outliers.tsv");
    std::string gene,species;
    std::vector<std::pair<std::string,std::string>> expected;
    while (pairs >> gene >> species) expected.emplace_back(gene,species);
    require(result.outliers==expected,"outliers or removal order differ from R");
    const auto cells = numbers(root/"carnivora.cells.tsv");
    require(cells.size()==2*result.cells.size(),"wrong cell count");
    for (std::size_t i = 0; i < result.cells.size(); ++i) {
      require(static_cast<double>(result.cells[i].first+1)==cells[2*i] &&
              static_cast<double>(result.cells[i].second+1)==cells[2*i+1],"removed cell order mismatch");
    }
    auto mc = input(root/"medcouple.tsv"); std::string line;
    while (std::getline(mc,line)) {
      std::istringstream row(line); double expected_mc,value; row >> expected_mc;
      std::vector<double> sample; while (row >> value) sample.push_back(value);
      require(std::abs(medcouple(sample)-expected_mc)<1e-13,"R/C++ medcouple mismatch");
    }
    for (std::size_t sample = 1; sample <= 10; ++sample) {
      const auto data = numbers(root/("cluster-"+std::to_string(sample)+".tsv"));
      Matrix matrix(12,12);
      require(data.size()==144,"wrong clustering fixture size");
      std::copy(data.begin(),data.end(),matrix.data().begin());
      const auto expected_order = numbers(root/("cluster-"+std::to_string(sample)+".order.tsv"));
      const auto order = cluster_order(matrix);
      require(order.size()==expected_order.size(),"clustering size mismatch");
      for (std::size_t i = 0; i < order.size(); ++i)
        require(static_cast<double>(order[i]+1)==expected_order[i],"R/C++ clustering tie mismatch");
    }
    std::cout << "PASS frozen R reference: 94 outliers, all scores, initial/final matrices, medcouple\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
