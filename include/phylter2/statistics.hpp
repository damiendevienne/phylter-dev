#pragma once
#include "phylter2/matrix.hpp"
#include <vector>
namespace phylter2 {
double medcouple(const std::vector<double>& values);
double quantile(std::vector<double> values, double probability);
double outlier_threshold(const std::vector<double>& values, double k);
std::vector<std::size_t> cluster_order(const Matrix& distances);
}
