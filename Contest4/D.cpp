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

class Graph {
 public:
  Graph() = default;

  Graph(long long vertex_numbers, long long edges_numbers)
      : edges_numbers_(edges_numbers), graph_(vertex_numbers + 1) {}

  void AddEdge(long long start, long long end, long long weight) {
    graph_[start].push_back({end, weight});
    graph_[end].push_back({start, weight});
  }

  const Vector1D<EdgeType>& GetNeighbors(long long vertex) const {
    return graph_[vertex];
  }

  long long GetVertexCount() const {
    return static_cast<long long>(graph_.size() - 1);
  }

  long long GetEdgesCount() const { return edges_numbers_; }

  const Vector2D<EdgeType>& GetAdjacencyList() const { return graph_; }

 private:
  long long edges_numbers_ = 0;
  Vector2D<EdgeType> graph_;
};

std::istream& operator>>(std::istream& stream, Graph& graph) {
  for (long long i = 0; i < graph.GetEdgesCount(); ++i) {
    long long start;
    long long end;
    long long weight;
    stream >> start >> end >> weight;
    graph.AddEdge(start, end, weight);
  }

  return stream;
}

class PrimsAlgorithm {
 public:
  PrimsAlgorithm() = default;

  PrimsAlgorithm(const Graph& graph) : graph_(graph) {}

  Graph CalculateMaxSpanningTree() {
    long long vertex_numbers = graph_.GetVertexCount();
    Graph max_spanning_tree(vertex_numbers, 0);

    std::set<EdgeType, std::greater<EdgeType>> vertex_queue;
    Vector1D<long long> max_edges(vertex_numbers + 1, kMinValue);
    Vector1D<long long> chosen_neighbors(vertex_numbers + 1, 0);
    std::vector<bool> uncaptures_vertices(vertex_numbers + 1, false);

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
        max_spanning_tree.AddEdge(best_neightbor, current_vertex,
                                  max_edges[current_vertex]);
      }

      for (const auto& [neighbor, weight] :
           graph_.GetNeighbors(current_vertex)) {
        if (!uncaptures_vertices[neighbor] && weight > max_edges[neighbor]) {
          vertex_queue.erase({max_edges[neighbor], neighbor});
          max_edges[neighbor] = weight;
          chosen_neighbors[neighbor] = current_vertex;
          vertex_queue.insert({weight, neighbor});
        }
      }
    }

    return max_spanning_tree;
  }

 private:
  Graph graph_;
};

class GetMinValueOnPath {
 public:
  GetMinValueOnPath() = default;

  GetMinValueOnPath(const Graph& graph) : graph_(graph) {
    long long vertex_numbers = graph_.GetVertexCount();
    log_vertex_numbers_ = std::ceil(std::log2(vertex_numbers)) + 1;
    depths_.resize(vertex_numbers);

    binary_ancestors_.resize(vertex_numbers,
                             Vector1D<long long>(log_vertex_numbers_, 0));
    min_binary_edges_.resize(
        vertex_numbers, Vector1D<long long>(log_vertex_numbers_, kMaxValue));

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

    for (const auto& [neighbor, weight] : graph_.GetNeighbors(vertex)) {
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

    for (const auto& [neighbor, weight] : graph_.GetNeighbors(vertex)) {
      if (neighbor != parent) {
        InitVertexInformation(neighbor, vertex, depth + 1);
      }
    }
  }

  Graph graph_;
  long long log_vertex_numbers_;

  Vector2D<long long> binary_ancestors_;
  Vector2D<long long> min_binary_edges_;
  Vector1D<long long> depths_;
};

int main() {
  long long vertex_numbers;
  long long edges_numbers;
  long long queries_numbers;
  std::cin >> vertex_numbers >> edges_numbers >> queries_numbers;

  Graph network(vertex_numbers, edges_numbers);
  std::cin >> network;

  PrimsAlgorithm prim(network);
  Graph max_spanning_tree_in_network = prim.CalculateMaxSpanningTree();

  GetMinValueOnPath calculator(max_spanning_tree_in_network);
  for (long long i = 0; i < queries_numbers; ++i) {
    long long server_start;
    long long server_end;
    std::cin >> server_start >> server_end;
    std::cout << calculator.GetMinEdgeValue(server_start, server_end) << "\n";
  }

  return 0;
}
