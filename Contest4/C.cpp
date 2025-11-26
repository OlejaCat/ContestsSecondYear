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

class PrimsAlgorithm {
 public:
  PrimsAlgorithm() = default;

  PrimsAlgorithm(std::size_t vertex_numbers)
      : vertex_numbers_(vertex_numbers),
        graph_(vertex_numbers + 1, Vector1D<EdgeType>()) {}

  void InitializeGraph(const Vector2D<std::size_t>& meeting_costs,
                       const Vector1D<std::size_t>& hiring_costs) {
    for (std::size_t i = 1; i <= vertex_numbers_; ++i) {
      for (std::size_t j = 1; j <= vertex_numbers_; ++j) {
        std::size_t weight = meeting_costs[i - 1][j - 1];
        graph_[i].push_back({j, weight});
        graph_[j].push_back({i, weight});
      }
    }

    for (std::size_t i = 1; i <= vertex_numbers_; ++i) {
      std::size_t weight = hiring_costs[i - 1];
      graph_[0].push_back({i, weight});
      graph_[i].push_back({0, weight});
    }
  }

  std::size_t CalculateMinSum() {
    std::set<EdgeType> vertex_queue;
    Vector1D<std::size_t> min_edges(vertex_numbers_ + 1,
                                    std::numeric_limits<std::size_t>::max());
    std::vector<bool> uncaptures_vertices(vertex_numbers_ + 1, false);

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

      for (const auto& [neighbor, weight] : graph_[current_vertex]) {
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
  std::size_t vertex_numbers_;
  Vector2D<EdgeType> graph_;
};

void FindCandidates();

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(0);

  FindCandidates();

  return 0;
}

void FindCandidates() {
  std::size_t candidates = 0;
  std::cin >> candidates;

  Vector2D<std::size_t> meeting_costs(candidates,
                                      Vector1D<std::size_t>(candidates));
  for (std::size_t i = 0; i < candidates; ++i) {
    for (std::size_t j = 0; j < candidates; ++j) {
      std::cin >> meeting_costs[i][j];
    }
  }

  Vector1D<std::size_t> hiring_costs(candidates);
  for (std::size_t i = 0; i < candidates; ++i) {
    std::cin >> hiring_costs[i];
  }

  PrimsAlgorithm prim(candidates);
  prim.InitializeGraph(meeting_costs, hiring_costs);
  std::cout << prim.CalculateMinSum() << "\n";
}
