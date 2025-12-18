#include <cstddef>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

template <typename T>
using Vector1D = std::vector<T>;
template <typename T>
using Vector2D = Vector1D<Vector1D<T>>;
using EdgeType = std::pair<std::size_t, std::size_t>;

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
  std::size_t vertex_count = graph.GetVertexCount();
  for (std::size_t i = 1; i <= vertex_count; ++i) {
    for (std::size_t j = 1; j <= vertex_count; ++j) {
      long long weight;
      stream >> weight;
      graph.AddEdge(i, j, weight);
    }
  }

  for (std::size_t i = 1; i <= vertex_count; ++i) {
    long long weight;
    stream >> weight;
    graph.AddEdge(0, i, weight);
  }

  return stream;
}

class PrimsAlgorithm {
 public:
  PrimsAlgorithm() = default;

  PrimsAlgorithm(const Graph& graph) : graph_(graph) {}

  std::size_t CalculateMinSum() {
    long long vertex_numbers = graph_.GetVertexCount();

    std::set<EdgeType> vertex_queue;
    Vector1D<std::size_t> min_edges(vertex_numbers + 1,
                                    std::numeric_limits<std::size_t>::max());
    std::vector<bool> uncaptures_vertices(vertex_numbers + 1, false);

    min_edges[0] = 0;
    vertex_queue.insert({0, 0});

    std::size_t min_weight = 0;

    while (!vertex_queue.empty()) {
      auto it = vertex_queue.begin();
      std::size_t current_vertex = it->second;
      std::size_t current_weight = it->first;
      vertex_queue.erase(it);

      if (uncaptures_vertices[current_vertex]) {
        continue;
      }

      uncaptures_vertices[current_vertex] = true;
      min_weight += current_weight;

      for (const auto& [neighbor, weight] :
           graph_.GetNeighbors(current_vertex)) {
        if (!uncaptures_vertices[neighbor] && weight < min_edges[neighbor]) {
          vertex_queue.erase({min_edges[neighbor], neighbor});
          min_edges[neighbor] = weight;
          vertex_queue.insert({weight, neighbor});
        }
      }
    }

    return min_weight;
  }

 private:
  Graph graph_;
};

std::size_t CalculateEdgesByCandidates(std::size_t candidates) {
  return candidates * (candidates - 1) + candidates;
}

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(0);

  std::size_t candidates = 0;
  std::cin >> candidates;

  Graph graph(candidates + 1, CalculateEdgesByCandidates(candidates));
  std::cin >> graph;

  PrimsAlgorithm prim(graph);
  std::cout << prim.CalculateMinSum() << "\n";

  return 0;
}
