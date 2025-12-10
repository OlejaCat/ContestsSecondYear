#include <algorithm>
#include <iostream>
#include <vector>

enum class Color { WHITE, GRAY, BLACK };

class DFSVisitor {
 public:
  virtual void Visit(int vertex) = 0;

  virtual void PostVisit(int vertex) = 0;
};

template <typename GraphT>
void RunDFS(int vertex, const GraphT& graph,
            std::vector<Color>& vertices_colors, DFSVisitor& visitor) {
  vertices_colors[vertex] = Color::GRAY;
  visitor.Visit(vertex);

  for (const auto& neighbor : graph[vertex]) {
    if (vertices_colors[neighbor] == Color::WHITE) {
      RunDFS(neighbor, graph, vertices_colors, visitor);
    }
  }

  vertices_colors[vertex] = Color::BLACK;
  visitor.PostVisit(vertex);
}

class ForwardPassVisitor : public DFSVisitor {
 public:
  ForwardPassVisitor(std::vector<int>& exit_order) : exit_order_(exit_order) {}

  void Visit(int vertex) override {}

  void PostVisit(int vertex) override { exit_order_.push_back(vertex); }

 private:
  std::vector<int>& exit_order_;
};

class ReversedPassVisitor : public DFSVisitor {
 public:
  ReversedPassVisitor(std::vector<int>& component)
      : current_component_(component) {}

  void Visit(int vertex) override { current_component_.push_back(vertex); }

  void PostVisit(int vertex) override {}

 private:
  std::vector<int>& current_component_;
};

class Graph {
 private:
  using GraphT = std::vector<std::vector<long long>>;

 public:
  Graph() = default;

  Graph(std::size_t vertices_number, std::size_t edges_number)
      : vertices_number_(vertices_number),
        edges_number_(edges_number),
        graph_(vertices_number + 1),
        graph_reversed_(vertices_number + 1) {}

  auto GetVerticesNumber() const { return vertices_number_; }
  auto GetEdgesNumber() const { return edges_number_; }

  const GraphT& GetForwardGraph() const { return graph_; }
  const GraphT& GetReversedGraph() const { return graph_reversed_; }

  void AddEdge(std::size_t lhs, std::size_t rhs) {
    graph_[lhs].push_back(rhs);
    graph_reversed_[rhs].push_back(lhs);
  }

 private:
  std::size_t vertices_number_ = 0;
  std::size_t edges_number_ = 0;

  GraphT graph_;
  GraphT graph_reversed_;
};

class StrongComponentsFinder {
 public:
  StrongComponentsFinder(const Graph& graph)
      : graph_(graph),
        vertices_colors_(graph.GetVerticesNumber() + 1, Color::WHITE),
        components_ids_(graph.GetVerticesNumber() + 1) {}

  void FindComponents() {
    FindStronglyConnectedComponents();  // Можно было бы и сюда реализацию, но
                                        // мне кажется так красивее
  }

  void PrintResult() const {
    std::cout << components_count_ << "\n";
    for (std::size_t i = 1; i <= graph_.GetVerticesNumber(); ++i) {
      std::cout << components_ids_[i] << " ";
    }
    std::cout << "\n";
  }

 private:
  void FindStronglyConnectedComponents() {
    int vertices_number = graph_.GetVerticesNumber();

    ForwardPassVisitor forward_visitor(exit_order_);
    for (int i = 1; i <= vertices_number; ++i) {
      if (vertices_colors_[i] == Color::WHITE) {
        RunDFS(i, graph_.GetForwardGraph(), vertices_colors_, forward_visitor);
      }
    }

    vertices_colors_.assign(vertices_number + 1, Color::WHITE);
    std::reverse(exit_order_.begin(), exit_order_.end());

    for (const auto& vertex : exit_order_) {
      if (vertices_colors_[vertex] == Color::WHITE) {
        current_component_.clear();

        ReversedPassVisitor reversed_visitor(current_component_);
        RunDFS(vertex, graph_.GetReversedGraph(), vertices_colors_,
               reversed_visitor);

        components_count_++;

        for (const auto& comp_vertex : current_component_) {
          components_ids_[comp_vertex] = components_count_;
        }
      }
    }
  }

  const Graph& graph_;

  std::vector<Color> vertices_colors_;
  std::vector<int> exit_order_;
  std::vector<int> current_component_;
  std::vector<int> components_ids_;

  int components_count_ = 0;
};

inline std::istream& operator>>(std::istream& istream, Graph& graph) {
  for (std::size_t i = 0; i < graph.GetEdgesNumber(); ++i) {
    std::size_t start;
    std::size_t end;
    std::cin >> start >> end;
    graph.AddEdge(start, end);
  }
  return istream;
}

int main() {
  std::size_t vertices_number;
  std::size_t edges_number;
  std::cin >> vertices_number >> edges_number;

  Graph graph(vertices_number, edges_number);
  std::cin >> graph;

  StrongComponentsFinder components_finder(graph);
  components_finder.FindComponents();

  components_finder.PrintResult();

  return 0;
}
