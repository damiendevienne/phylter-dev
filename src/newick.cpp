#include "phylter2/newick.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace phylter2 {
namespace {

class Parser {
 public:
  explicit Parser(const std::string_view source) : source_(source) {}

  Tree parse() {
    auto result = parse_one();
    if (position_ != source_.size()) throw error("unexpected text after ';'");
    return result;
  }

  std::vector<Tree> parse_many() {
    std::vector<Tree> result;
    skip_ignored();
    while (position_ < source_.size()) result.push_back(parse_one());
    if (result.empty()) throw error("no trees found");
    return result;
  }

 private:
  Tree parse_one() {
    tree_ = Tree();
    skip_ignored();
    const std::size_t root = subtree();
    skip_ignored();
    expect(';');
    skip_ignored();
    tree_.root = root;
    return std::move(tree_);
  }

 private:
  [[nodiscard]] std::runtime_error error(const std::string& message) const {
    return std::runtime_error("invalid Newick at byte " + std::to_string(position_ + 1) + ": " + message);
  }

  void skip_ignored() {
    while (position_ < source_.size()) {
      if (std::isspace(static_cast<unsigned char>(source_[position_]))) {
        ++position_;
        continue;
      }
      if (source_[position_] != '[') return;
      ++position_;
      int nesting = 1;
      while (position_ < source_.size() && nesting > 0) {
        if (source_[position_] == '[') ++nesting;
        if (source_[position_] == ']') --nesting;
        ++position_;
      }
      if (nesting != 0) throw error("unterminated bracket comment");
    }
  }

  void expect(const char expected) {
    skip_ignored();
    if (position_ == source_.size() || source_[position_] != expected) {
      throw error(std::string("expected '") + expected + "'");
    }
    ++position_;
  }

  [[nodiscard]] bool next_is(const char expected) {
    skip_ignored();
    return position_ < source_.size() && source_[position_] == expected;
  }

  [[nodiscard]] std::string label() {
    skip_ignored();
    if (position_ == source_.size()) return {};
    if (source_[position_] == '\'') {
      ++position_;
      std::string value;
      while (position_ < source_.size()) {
        const char character = source_[position_++];
        if (character != '\'') {
          value += character;
        } else if (position_ < source_.size() && source_[position_] == '\'') {
          value += '\'';
          ++position_;
        } else {
          return value;
        }
      }
      throw error("unterminated quoted label");
    }
    const std::size_t begin = position_;
    while (position_ < source_.size()) {
      const char character = source_[position_];
      if (std::isspace(static_cast<unsigned char>(character)) || character == ':' || character == ',' ||
          character == '(' || character == ')' || character == ';' || character == '[') {
        break;
      }
      ++position_;
    }
    return std::string(source_.substr(begin, position_ - begin));
  }

  void optional_label_and_length(TreeNode& node) {
    node.label = label();
    skip_ignored();
    if (!next_is(':')) return;
    ++position_;
    skip_ignored();
    // strtod requires a terminated buffer; string_view does not promise one.
    auto finish = position_;
    while (finish < source_.size() && source_[finish] != ',' && source_[finish] != ')' &&
           source_[finish] != ';' && source_[finish] != '[' &&
           !std::isspace(static_cast<unsigned char>(source_[finish]))) ++finish;
    const std::string remaining(source_.substr(position_,finish-position_));
    const char* begin = remaining.c_str();
    char* end = nullptr;
    const double value = std::strtod(begin, &end);
    if (end == begin || !std::isfinite(value) || value < 0.0) {
      throw error("edge length must be a finite non-negative number");
    }
    position_ += static_cast<std::size_t>(end - begin);
    node.length_to_parent = value;
    node.has_length = true;
  }

  [[nodiscard]] std::size_t subtree() {
    if (++depth_ > 2048) throw error("tree nesting exceeds the supported depth (2048)");
    skip_ignored();
    const std::size_t id = tree_.nodes.size();
    tree_.nodes.emplace_back();
    if (next_is('(')) {
      ++position_;
      const std::size_t first_child = subtree();
      tree_.nodes[id].children.push_back(first_child);
      while (next_is(',')) {
        ++position_;
        const std::size_t child = subtree();
        tree_.nodes[id].children.push_back(child);
      }
      expect(')');
      optional_label_and_length(tree_.nodes[id]);
      --depth_;
      return id;
    }
    optional_label_and_length(tree_.nodes[id]);
    if (tree_.nodes[id].label.empty()) throw error("leaf label is required");
    --depth_;
    return id;
  }

  std::string_view source_;
  std::size_t position_{0};
  std::size_t depth_{0};
  Tree tree_;
};

}  // namespace

Tree parse_newick(const std::string_view source) {
  return Parser(source).parse();
}

std::vector<Tree> parse_newick_set(const std::string_view source) {
  return Parser(source).parse_many();
}

NamedDistanceMatrix patristic_distances(const Tree& tree, const bool nodal, const double support_cutoff) {
  if (!std::isfinite(support_cutoff) || support_cutoff < 0)
    throw std::invalid_argument("support cutoff must be finite and non-negative");
  if (tree.nodes.empty() || tree.root >= tree.nodes.size()) {
    throw std::invalid_argument("tree has no valid root");
  }
  std::vector<std::vector<std::pair<std::size_t, double>>> edges(tree.nodes.size());
  std::vector<std::size_t> parents(tree.nodes.size());
  std::vector<std::size_t> leaves;
  std::vector<std::string> taxa;
  std::unordered_set<std::string> observed;
  const bool has_support = std::any_of(tree.nodes.begin(),tree.nodes.end(),[](const auto& node) {
    return !node.children.empty() && !node.label.empty();
  });
  for (std::size_t parent = 0; parent < tree.nodes.size(); ++parent) {
    const TreeNode& node = tree.nodes[parent];
    bool low_support = false;
    if (support_cutoff != 0 && has_support && !node.children.empty()) {
      std::size_t end = 0;
      const double support = std::stod(node.label,&end);
      if (end != node.label.size() || !std::isfinite(support))
        throw std::invalid_argument("support filtering requires numeric labels on all internal nodes");
      low_support = support < support_cutoff;
    }
    if (node.children.empty()) {
      if (node.label.empty() || !observed.insert(node.label).second) {
        throw std::invalid_argument("tree leaves must have unique non-empty labels");
      }
      leaves.push_back(parent);
      taxa.push_back(node.label);
    }
    for (const std::size_t child : node.children) {
      if (child >= tree.nodes.size() || (!nodal && !tree.nodes[child].has_length)) {
        throw std::invalid_argument("every non-root tree edge needs a branch length");
      }
      if (child == tree.root || ++parents[child] != 1)
        throw std::invalid_argument("tree nodes must have exactly one parent");
      double length = tree.nodes[child].length_to_parent;
      const bool collapse = support_cutoff != 0 && has_support;
      if (low_support) length = 1e-10;
      // Preserve trees2matrices: low-support PARENTS modify outgoing edges;
      // di2multi collapses all short INTERNAL edges, never terminal edges.
      const bool collapsed = collapse && !tree.nodes[child].children.empty() && length < 1e-9;
      if (nodal) length = collapsed ? 0.0 : 1.0;
      else if (collapsed) length = 0.0;
      if (!std::isfinite(length) || length < 0)
        throw std::invalid_argument("invalid edge length");
      edges[parent].emplace_back(child, length);
      edges[child].emplace_back(parent, length);
    }
  }
  for (std::size_t i = 0; i < parents.size(); ++i)
    if (i != tree.root && parents[i] != 1) throw std::invalid_argument("disconnected tree node");
  std::vector<std::size_t> pending{tree.root};
  std::size_t reachable = 0;
  while (!pending.empty()) {
    const auto node = pending.back(); pending.pop_back(); ++reachable;
    for (auto child : tree.nodes[node].children) pending.push_back(child);
  }
  if (reachable != tree.nodes.size()) throw std::invalid_argument("disconnected cycle in tree");
  if (leaves.size() < 3) throw std::invalid_argument("tree requires at least three labelled leaves");

  Matrix distances(leaves.size(), leaves.size());
  for (std::size_t source = 0; source < leaves.size(); ++source) {
    std::vector<double> distance(tree.nodes.size(), std::numeric_limits<double>::infinity());
    std::vector<std::size_t> stack{leaves[source]};
    distance[leaves[source]] = 0.0;
    while (!stack.empty()) {
      const std::size_t node = stack.back();
      stack.pop_back();
      for (const auto [neighbor, length] : edges[node]) {
        if (std::isfinite(distance[neighbor])) continue;
        distance[neighbor] = distance[node] + length;
        stack.push_back(neighbor);
      }
    }
    for (std::size_t target = source + 1; target < leaves.size(); ++target) {
      if (!std::isfinite(distance[leaves[target]])) throw std::invalid_argument("disconnected tree");
      distances(source, target) = distances(target, source) = distance[leaves[target]] - (nodal ? 1.0 : 0.0);
    }
  }
  return {.values = std::move(distances), .taxa = std::move(taxa)};
}

}  // namespace phylter2
