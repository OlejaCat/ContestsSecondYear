#include <algorithm>
#include <iostream>
#include <vector>

class Solution {
 private:
  enum class Color { WHITE, GRAY, BLACK };

 public:
  Solution() = default;

  void Solve() {
    int vertices_number = 0;
    int edges_number = 0;
    std::cin >> vertices_number >> edges_number;

    graph_.resize(vertices_number + 1);
    graph_reversed_.resize(vertices_number + 1);
    vertices_colors_.assign(vertices_number + 1, Color::WHITE);
    components_ids_.resize(vertices_number + 1);

    BuildGraph(edges_number);
    FindStronglyConnectedComponents(vertices_number);

    PrintResult(vertices_number);
  }

 private:
  void PrintResult(int vertices_number) {
    std::cout << components_count_ << "\n";
    for (int i = 1; i <= vertices_number; ++i) {
      std::cout << components_ids_[i] << " ";
    }
    std::cout << "\n";
  }

  void BuildGraph(int edges_number) {
    for (int i = 0; i < edges_number; ++i) {
      int start = 0;
      int end = 0;
      std::cin >> start >> end;
      graph_[start].push_back(end);
      graph_reversed_[end].push_back(start);
    }
  }

  void FindStronglyConnectedComponents(int vertices_number) {
    for (int i = 1; i <= vertices_number; ++i) {
      if (vertices_colors_[i] == Color::WHITE) {
        ForwardDFS(i);
      }
    }

    vertices_colors_.assign(vertices_number + 1, Color::WHITE);
    std::reverse(exit_order_.begin(), exit_order_.end());

    for (const auto& vertex : exit_order_) {
      if (vertices_colors_[vertex] == Color::WHITE) {
        current_component_.clear();
        ReversedDFS(vertex);
        components_count_++;

        for (const auto& comp_vertex : current_component_) {
          components_ids_[comp_vertex] = components_count_;
        }
      }
    }
  }

  void ForwardDFS(int vertex) {
    vertices_colors_[vertex] = Color::GRAY;

    for (const auto& adj_vertex : graph_[vertex]) {
      if (vertices_colors_[adj_vertex] == Color::WHITE) {
        ForwardDFS(adj_vertex);
      }
    }

    vertices_colors_[vertex] = Color::BLACK;
    exit_order_.push_back(vertex);
  }

  void ReversedDFS(int vertex) {
    vertices_colors_[vertex] = Color::GRAY;
    current_component_.push_back(vertex);

    for (const auto& adj_vertex : graph_reversed_[vertex]) {
      if (vertices_colors_[adj_vertex] == Color::WHITE) {
        ReversedDFS(adj_vertex);
      }
    }

    vertices_colors_[vertex] = Color::BLACK;
  }

  std::vector<std::vector<int>> graph_;
  std::vector<std::vector<int>> graph_reversed_;
  std::vector<Color> vertices_colors_;
  std::vector<int> exit_order_;
  std::vector<int> current_component_;
  std::vector<int> components_ids_;

  int components_count_ = 0;
};

int main() {
  Solution solution;
  solution.Solve();

  return 0;
}
