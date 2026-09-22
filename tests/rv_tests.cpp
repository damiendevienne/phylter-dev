#include "phylter2/rv.hpp"
#include "phylter2/lapack.hpp"
#include "phylter2/distatis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace phylter2;
namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
std::vector<double> multiply(const Matrix& a,const std::vector<double>& x) {
  std::vector<double> y(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    for (std::size_t j = 0; j < x.size(); ++j) y[i] += a(i,j)*x[j];
  return y;
}
void near(double a,double b,double tolerance,const char* message) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=tolerance,message);
}
template<class Function> void rejects(Function function) {
  try { function(); } catch(const std::exception&) { return; }
  throw std::runtime_error("invalid or nonconverged eigenproblem accepted");
}
}
int main() {
  try {
    // The implicit operator must reproduce the upper-triangle convention even
    // when lower entries differ (legacy nearest-neighbor imputation).
    std::vector<Matrix> matrices;
    for (std::size_t gene = 0; gene < 70; ++gene) {
      Matrix matrix(5,5);
      for (std::size_t i = 0; i < 5; ++i) for (std::size_t j = 0; j < 5; ++j)
        matrix(i,j) = std::sin(static_cast<double>(gene*17+i*7+j*11+1));
      matrices.push_back(std::move(matrix));
    }
    RvOperator op(matrices);
    const auto dense = op.materialize();
    for (std::size_t sample = 0; sample < 4; ++sample) {
      std::vector<double> x(matrices.size());
      for (std::size_t i = 0; i < x.size(); ++i) x[i] = std::cos(static_cast<double>(i*7+sample));
      const auto actual = op(x), expected = multiply(dense,x);
      for (std::size_t i = 0; i < x.size(); ++i)
        near(actual[i],expected[i],2e-13,"implicit RV product differs from dense RV");
    }
    const auto iterative = leading_eigenpair(matrices.size(),[&](const auto& x){return op(x);});
    const auto exact = symmetric_eigen(dense);
    near(iterative.value,exact.values.front(),1e-11,"RV eigenvalue mismatch");
    require(iterative.relative_residual<=2e-14,"RV residual exceeds tolerance");
    const auto image = multiply(dense,iterative.vector);
    for (std::size_t i = 0; i < image.size(); ++i)
      near(image[i],iterative.value*iterative.vector[i],1e-11,"dense RV residual check failed");

    // Rank-one RV, with a large null eigenspace.
    std::vector<Matrix> duplicates(130,matrices[0]);
    RvOperator duplicate_op(duplicates);
    const auto repeated = leading_eigenpair(duplicates.size(),[&](const auto& x){return duplicate_op(x);});
    near(repeated.value,130,1e-11,"rank-one RV failed");

    // Constant initialization is EXACTLY orthogonal to the dominant vector.
    // A plain power iteration from ones would falsely return eigenvalue 1.
    Matrix orthogonal(64,64);
    for (std::size_t i = 0; i < 64; ++i) for (std::size_t j = 0; j < 64; ++j)
      orthogonal(i,j) = (i==j ? 1.0 : 0.0) + 3.0*(i%2==j%2 ? 1.0 : -1.0)/64;
    const auto hidden = leading_eigenpair(64,[&](const auto& x){return multiply(orthogonal,x);});
    near(hidden.value,4,1e-12,"solver missed dominant vector orthogonal to constant start");

    // Force multiple restarts with a small subspace.
    Matrix diagonal(100,100);
    for (std::size_t i = 0; i < 100; ++i) diagonal(i,i) = .01*static_cast<double>(i);
    diagonal(99,99)=1.1;
    const auto restarted = leading_eigenpair(100,[&](const auto& x){return multiply(diagonal,x);},
                                            {.tolerance=2e-14,.max_basis=12,.max_products=2048});
    near(restarted.value,1.1,1e-12,"restarted solve failed");
    require(restarted.products>12,"restart test did not exercise restarts");
    diagonal(98,98)=1.1-1e-7;
    const auto clustered = leading_eigenpair(100,[&](const auto& x){return multiply(diagonal,x);});
    near(clustered.value,1.1,1e-12,"clustered leading eigenvalues failed");
    require(std::abs(clustered.vector[99])>0.999999,"wrong vector for clustered spectrum");
    diagonal(98,98)=1.1;
    const auto tied = leading_eigenpair(100,[&](const auto& x){return multiply(diagonal,x);});
    near(tied.value,1.1,1e-12,"tied leading eigenspace failed");

    rejects([&]{(void)leading_eigenpair(100,[&](const auto& x){return multiply(diagonal,x);},
                                      {.tolerance=2e-14,.max_basis=12,.max_products=1});});
    rejects([]{(void)leading_eigenpair(2,[](const auto&){return std::vector<double>{0};});});
    rejects([]{(void)leading_eigenpair(2,[](const auto&){return std::vector<double>{0,std::numeric_limits<double>::quiet_NaN()};});});
    rejects([]{std::vector<Matrix> zero{Matrix(3,3)}; (void)RvOperator(zero);});

    // Same full DISTATIS problem with and without the stored RV.
    std::vector<Matrix> distances;
    for (std::size_t gene = 0; gene < 9; ++gene) {
      Matrix d(8,8);
      for (std::size_t i = 0; i < 8; ++i) for (std::size_t j = i+1; j < 8; ++j)
        d(i,j)=d(j,i)=1+std::abs(std::sin(static_cast<double>(i*13+j*3+gene*7)));
      distances.push_back(std::move(d));
    }
    const auto baseline = distatis(distances,0,-1,{RvMethod::dense,true});
    const auto implicit = distatis(distances,0,-1,{RvMethod::matrix_free,false});
    require(implicit.rv.rows()==0 && implicit.rv.cols()==0,"matrix-free mode retained RV");
    near(implicit.quality,baseline.quality,1e-12,"DISTATIS quality changed");
    require(max_abs_difference(dist2wr(implicit),dist2wr(baseline))<1e-10,"DISTATIS WR changed");
    require(max_abs_difference(implicit.compromise,baseline.compromise)<1e-12,"compromise changed");
    const auto rejected = distatis(distances,0,2,{RvMethod::matrix_free,false});
    require(rejected.rv.rows()==0 && rejected.partial_coordinates.empty(),"rejected proposal retained matrices");
    // Large gene count with few taxa: storing/eigendecomposing a 5000x5000 RV
    // would make this small regression unacceptably expensive.
    std::vector<Matrix> many(5000,distances.front());
    const auto scalable = distatis(many,0,-1,{RvMethod::matrix_free,false});
    near(scalable.quality,1,1e-11,"many-gene rank-one quality changed");
    require(scalable.rv.rows()==0,"many-gene solve materialized RV");
    std::cout << "PASS implicit RV, dense equivalence, restarted eigensolver and difficult spectra\n";
  } catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
