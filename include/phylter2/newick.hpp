#pragma once

#include "phylter2/matrix.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace phylter2 {

struct TreeNode {
  std::string label;
  std::vector<std::size_t> children;
  double length_to_parent{0.0};
  bool has_length{false};
};

struct Tree {
  std::vector<TreeNode> nodes;
  std::size_t root{0};
};

struct NamedDistanceMatrix {
  Matrix values;
  std::vector<std::string> taxa;
};

// Parses a standard Newick tree. Quoted labels and bracket comments are
// accepted. A string may contain exactly one tree, terminated by `;`.
[[nodiscard]] Tree parse_newick(std::string_view source);
[[nodiscard]] std::vector<Tree> parse_newick_set(std::string_view source);

// Computes pairwise sums of branch lengths between leaves. All non-root edges
// must have finite, non-negative lengths and all leaves need unique labels.
[[nodiscard]] NamedDistanceMatrix patristic_distances(const Tree& tree, bool nodal = false,
                                                      double support_cutoff = 0);

}  // namespace phylter2
