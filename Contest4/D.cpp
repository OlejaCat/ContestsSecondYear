#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

template <typename T>
using Vector1D = std::vector<T>;
template <typename T>
using Vector2D = Vector1D<Vector1D<T>>;
using EdgeType = std::pair<long long, long long>;

const long long kMinValue = std::numeric_limits<long long>::min();
const long long kMaxValue = std::numeric_limits<long long>::max();

class PrimsAlgorithm {
 public:
  PrimsAlgorithm() = default;

  PrimsAlgorithm(long long vertex_numbers, long long edges_numbers)
      : vertex_numbers_(vertex_numbers),
        edges_numbers_(edges_numbers),
        graph_(vertex_numbers + 1, Vector1D<EdgeType>()) {}

  void InitializeGraph() {
    for (long long edge = 0; edge < edges_numbers_; ++edge) {
      long long start;
      long long end;
      long long weight;
      std::cin >> start >> end >> weight;
      graph_[start].push_back({end, weight});
      graph_[end].push_back({start, weight});
    }
  }

  void CalculateMaxSpanningTree() {
    std::set<EdgeType, std::greater<EdgeType>> vertex_queue;
    Vector1D<long long> max_edges(vertex_numbers_ + 1, kMinValue);
    Vector1D<long long> chosen_neighbors(vertex_numbers_ + 1, 0);
    std::vector<bool> uncaptures_vertices(vertex_numbers_ + 1, false);

    max_edges[1] = kMinValue;
    vertex_queue.insert({kMinValue, 1});

    while (!vertex_queue.empty()) {
      auto it = vertex_queue.begin();
      long long current_vertex = it->second;
      vertex_queue.erase(it);

      if (uncaptures_vertices[current_vertex]) {
        continue;
      }

      uncaptures_vertices[current_vertex] = true;

      long long best_neightbor = chosen_neighbors[current_vertex];
      if (best_neightbor != 0) {
        max_spanning_tree_[best_neightbor].push_back(
            {current_vertex, max_edges[current_vertex]});
        max_spanning_tree_[current_vertex].push_back(
            {best_neightbor, max_edges[current_vertex]});
      }

      for (const auto& [neighbor, weight] : graph_[current_vertex]) {
        if (!uncaptures_vertices[neighbor] && weight > max_edges[neighbor]) {
          vertex_queue.erase({max_edges[neighbor], neighbor});
          max_edges[neighbor] = weight;
          chosen_neighbors[neighbor] = current_vertex;
          vertex_queue.insert({weight, neighbor});
        }
      }
    }
  }

  Vector2D<EdgeType>& GetTree() {
    max_spanning_tree_.clear();
    max_spanning_tree_.assign(vertex_numbers_ + 1, Vector1D<EdgeType>());
    CalculateMaxSpanningTree();
    return max_spanning_tree_;
  }

 private:
  long long vertex_numbers_ = 0;
  long long edges_numbers_ = 0;
  Vector2D<EdgeType> graph_;
  Vector2D<EdgeType> max_spanning_tree_;
};

class MaxSpanningTree {
 public:
  MaxSpanningTree() = default;

  MaxSpanningTree(const Vector2D<EdgeType>& graph) : graph_(graph) {
    log_vertex_numbers_ = std::ceil(std::log2(graph_.size())) + 1;
    depths_.resize(graph_.size());

    binary_ancestors_.resize(graph_.size(),
                             Vector1D<long long>(log_vertex_numbers_, 0));
    min_binary_edges_.resize(
        graph_.size(), Vector1D<long long>(log_vertex_numbers_, kMaxValue));

    InitVertexInformation(1);
  }

  long long GetMinEdgeValue(long long lhs, long long rhs) {
    long long min_edge = kMaxValue;

    if (depths_[lhs] < depths_[rhs]) {
      std::swap(lhs, rhs);
    }

    for (ssize_t i = log_vertex_numbers_ - 1; i >= 0; --i) {
      if (depths_[lhs] - (1 << i) >= depths_[rhs]) {
        min_edge = std::min(min_edge, min_binary_edges_[lhs][i]);
        lhs = binary_ancestors_[lhs][i];
      }
    }

    if (lhs == rhs) {
      return min_edge;
    }

    for (ssize_t i = log_vertex_numbers_ - 1; i >= 0; --i) {
      if (binary_ancestors_[lhs][i] != binary_ancestors_[rhs][i]) {
        min_edge = std::min(min_edge, min_binary_edges_[lhs][i]);
        min_edge = std::min(min_edge, min_binary_edges_[rhs][i]);
        lhs = binary_ancestors_[lhs][i];
        rhs = binary_ancestors_[rhs][i];
      }
    }

    return std::min(
        {min_edge, min_binary_edges_[lhs][0], min_binary_edges_[rhs][0]});
  }

 private:
  long long GetWeightToParent(long long vertex, long long parent) {
    if (parent == 0) {
      return kMaxValue;
    }

    for (const auto& [neighbor, weight] : graph_[vertex]) {
      if (neighbor == parent) {
        return weight;
      }
    }

    return kMaxValue;
  }

  void InitVertexInformation(long long vertex, long long parent = 0,
                             long long depth = 0) {
    depths_[vertex] = depth;
    binary_ancestors_[vertex][0] = parent;
    min_binary_edges_[vertex][0] = GetWeightToParent(vertex, parent);

    for (long long i = 1; i < log_vertex_numbers_; ++i) {
      long long prev_ancestor = binary_ancestors_[vertex][i - 1];
      binary_ancestors_[vertex][i] = binary_ancestors_[prev_ancestor][i - 1];
      min_binary_edges_[vertex][i] =
          std::min(min_binary_edges_[vertex][i - 1],
                   min_binary_edges_[prev_ancestor][i - 1]);
    }

    for (const auto& [neighbor, weight] : graph_[vertex]) {
      if (neighbor != parent) {
        InitVertexInformation(neighbor, vertex, depth + 1);
      }
    }
  }

  Vector2D<EdgeType> graph_;
  long long log_vertex_numbers_;

  Vector2D<long long> binary_ancestors_;
  Vector2D<long long> min_binary_edges_;
  Vector1D<long long> depths_;
};

void AnswerQueries();

int main() {
  AnswerQueries();

  return 0;
}

void AnswerQueries() {
  long long vertex_numbers;
  long long edges_numbers;
  long long queries_numbers;
  std::cin >> vertex_numbers >> edges_numbers >> queries_numbers;

  PrimsAlgorithm prim(vertex_numbers, edges_numbers);
  prim.InitializeGraph();

  MaxSpanningTree tree(prim.GetTree());
  for (long long i = 0; i < queries_numbers; ++i) {
    long long lhs;
    long long rhs;
    std::cin >> lhs >> rhs;
    std::cout << tree.GetMinEdgeValue(lhs, rhs) << "\n";
  }
}
