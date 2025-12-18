#include <cassert>
#include <iostream>
#include <limits>
#include <vector>
#include <queue>
#include <cmath>
#include <set>
#include <memory>

const long long kMaxValue = std::numeric_limits<long long>::max();

template <typename T>
using Vector1D = std::vector<T>;

template <typename T>
using Vector2D = Vector1D<Vector1D<T>>;

class Graph {
 private:
  const long long kZeroDistance = 0;

 public:
  Graph() = default;

  Graph(long long vertex_numbers, long long edge_numbers)
      : edge_numbers_(edge_numbers),
        graph_(vertex_numbers + 1) {}

  Graph(const Graph& other)
      : edge_numbers_(other.edge_numbers_),
        graph_(other.graph_),
        distances_from_vertex_(other.distances_from_vertex_) {}

  Graph& operator=(const Graph& other) {
    if (this != &other) {
      edge_numbers_ = other.edge_numbers_;
      graph_ = other.graph_;
      distances_from_vertex_ = other.distances_from_vertex_;
    }
    return *this;
  }

  void AddEdge(long long start, long long end) {
    assert(start <= GetVertexNumber());
    assert(end <= GetVertexNumber());

    graph_[start].push_back(end);
    graph_[end].push_back(start);
  }

  Vector1D<long long>& GetNeighbors(long long vertex) {
    assert(vertex <= GetVertexNumber());

    return graph_[vertex];
  }

  const Vector1D<long long>& GetNeighbors(long long vertex) const {
    assert(vertex <= GetVertexNumber());

    return graph_[vertex];
  }

  long long GetVertexNumber() const { return graph_.size() - 1; }

  void CalculateDistances(long long finish_vertex) {
    assert(finish_vertex <= GetVertexNumber());

    distances_from_vertex_.assign(graph_.size(), kMaxValue);

    std::queue<long long> vertex_to_process;
    distances_from_vertex_[finish_vertex] = kZeroDistance;
    vertex_to_process.push(finish_vertex);

    while (!vertex_to_process.empty()) {
      long long current_vertex = vertex_to_process.front();
      vertex_to_process.pop();

      for (long long neighbor : GetNeighbors(current_vertex)) {
        if (distances_from_vertex_[neighbor] == kMaxValue) {
          distances_from_vertex_[neighbor] =
              distances_from_vertex_[current_vertex] + 1;
          vertex_to_process.push(neighbor);
        }
      }
    }
  }

  long long GetDistance(long long vertex) const {
    assert(vertex <= GetVertexNumber());

    return distances_from_vertex_[vertex];
  }

 private:
  long long edge_numbers_;
  Vector2D<long long> graph_;
  Vector1D<long long> distances_from_vertex_;
};


class BridgeFinder {
 private:
  const long long kStartingEdge = 1;

 public:
  BridgeFinder() = default;

  BridgeFinder(const Graph& graph) : graph_(graph) {}

  BridgeFinder(const BridgeFinder& other)
      : graph_(other.graph_),
        bridges_(other.bridges_),
        used_vertex_(other.used_vertex_),
        min_depths_inverse_(other.min_depths_inverse_),
        depths_(other.depths_) {}

  BridgeFinder& operator=(const BridgeFinder& other) {
    if (this != &other) {
      graph_ = other.graph_;
      bridges_ = other.bridges_;
      used_vertex_ = other.used_vertex_;
      min_depths_inverse_ = other.min_depths_inverse_;
      depths_ = other.depths_;
    }
    return *this;
  }

  void FindBridges() {
    long long vertex_numbers = graph_.GetVertexNumber();

    used_vertex_.assign(vertex_numbers + 1, false);
    min_depths_inverse_.assign(vertex_numbers + 1, 0);
    depths_.assign(vertex_numbers + 1, 0);

    for (long long i = 1; i <= vertex_numbers; ++i) {
      if (!used_vertex_[i]) {
        FindBridgesDFS(i);
      }
    }
  }

  bool IsBridge(long long lhs, long long rhs) const {
    return bridges_.count({std::min(lhs, rhs), std::max(lhs, rhs)}) > 0;
  }

 private:
  void FindBridgesDFS(long long vertex, long long parent = 0) {
    used_vertex_[vertex] = true;
    depths_[vertex] = (parent == 0) ? 0 : depths_[parent] + 1;
    min_depths_inverse_[vertex] = depths_[vertex];

    for (const auto& neighbor : graph_.GetNeighbors(vertex)) {
      if (neighbor == parent)
        continue;

      if (used_vertex_[neighbor]) {
        min_depths_inverse_[vertex] =
            std::min(min_depths_inverse_[vertex], depths_[neighbor]);
        continue;
      }

      FindBridgesDFS(neighbor, vertex);
      min_depths_inverse_[vertex] =
          std::min(min_depths_inverse_[vertex], min_depths_inverse_[neighbor]);

      if (min_depths_inverse_[neighbor] > depths_[vertex]) {
        bridges_.insert(
            {std::min(vertex, neighbor), std::max(vertex, neighbor)});
      }
    }
  }

  Graph graph_;
  std::set<std::pair<long long, long long>> bridges_;
  Vector1D<bool> used_vertex_;
  Vector1D<long long> min_depths_inverse_;
  Vector1D<long long> depths_;
};

class ComponentFinder {
 public:
  ComponentFinder(const Graph& graph, const BridgeFinder& bridge_finder)
      : graph_(graph), bridge_finder_(bridge_finder) {}

  void FindComponents() {
    long long vertex_numbers = graph_.GetVertexNumber();
    component_.assign(vertex_numbers + 1, 0);
    component_id_ = 0;

    for (long long i = 1; i <= vertex_numbers; ++i) {
      if (component_[i] == 0) {
        ++component_id_;
        FindComponentsDFS(i, component_id_);
      }
    }
  }

  long long GetComponentId(long long vertex) const {
    return component_[vertex];
  }

  long long GetComponentCount() const { return component_id_; }

 private:
  void FindComponentsDFS(long long vertex, long long id) {
    component_[vertex] = id;
    for (const auto& neighbor : graph_.GetNeighbors(vertex)) {
      if (component_[neighbor] == 0
          && !bridge_finder_.IsBridge(vertex, neighbor)) {
        FindComponentsDFS(neighbor, id);
      }
    }
  }

  Graph graph_;
  const BridgeFinder& bridge_finder_;
  Vector1D<long long> component_;
  long long component_id_;
};

class BridgeTree {
 public:
  BridgeTree(const Graph& graph, const ComponentFinder& component_finder)
      : graph_(graph), component_finder_(component_finder) {}

  BridgeTree(const BridgeTree& other) = delete;
  BridgeTree& operator=(const BridgeTree& other) = delete;

  void BuildBridgeTree(long long finish_vertex) {
    BuildTreeGraph();
    CalculateDistancesInBridgeTree(finish_vertex);
  }

  const Graph& GetTreeGraph() const { return bridge_tree_graph_; }

  const Vector1D<long long>& GetDistancesInBridgeTree() const {
    return distances_in_bridge_tree_;
  }

 private:
  void BuildTreeGraph() {
    long long components_number = component_finder_.GetComponentCount();
    bridge_tree_graph_ = Graph(components_number, 0);

    for (long long vertex = 1; vertex <= graph_.GetVertexNumber(); ++vertex) {
      for (const auto& neighbor : graph_.GetNeighbors(vertex)) {
        if (vertex < neighbor) {
          long long component_vertex_id =
              component_finder_.GetComponentId(vertex);
          long long component_neighbor_id =
              component_finder_.GetComponentId(neighbor);
          if (component_vertex_id != component_neighbor_id) {
            bridge_tree_graph_.AddEdge(component_vertex_id,
                                       component_neighbor_id);
          }
        }
      }
    }
  }

  void CalculateDistancesInBridgeTree(long long finish_vertex) {
    long long components_number = component_finder_.GetComponentCount();
    distances_in_bridge_tree_.assign(components_number + 1, kMaxValue);

    long long finish_component =
        component_finder_.GetComponentId(finish_vertex);
    std::queue<long long> queue;
    distances_in_bridge_tree_[finish_component] = 0;
    queue.push(finish_component);

    while (!queue.empty()) {
      long long current_component = queue.front();
      queue.pop();

      for (long long neighbor :
           bridge_tree_graph_.GetNeighbors(current_component)) {
        if (distances_in_bridge_tree_[neighbor] == kMaxValue) {
          distances_in_bridge_tree_[neighbor] =
              distances_in_bridge_tree_[current_component] + 1;
          queue.push(neighbor);
        }
      }
    }
  }

  const Graph& graph_;
  const ComponentFinder& component_finder_;
  Graph bridge_tree_graph_;
  Vector1D<long long> distances_in_bridge_tree_;
};

class LcaFinder {
 public:
  LcaFinder(const Graph& tree, long long root,
            const Vector1D<long long>& distances_in_bridge_tree)
      : tree_(tree),
        root_(root),
        distances_in_bridge_tree_(distances_in_bridge_tree) {
    Precompute();
  }

  LcaFinder(const LcaFinder& other) = delete;
  LcaFinder& operator=(const LcaFinder& other) = delete;

  void Precompute() {
    long long tree_size = tree_.GetVertexNumber();

    log_tree_size_ = 0;
    if (tree_size > 0) {
      log_tree_size_ =
          static_cast<long long>(std::ceil(std::log2(tree_size))) + 1;
    }

    depths_.assign(tree_size + 1, 0);
    binary_ancestors_.assign(tree_size + 1,
                             Vector1D<long long>(log_tree_size_, 0));

    DfsForLCA(root_, 0, 0);
  }

  long long FindLCA(long long lhs, long long rhs) const {
    if (depths_[lhs] < depths_[rhs]) {
      std::swap(lhs, rhs);
    }

    for (int i = log_tree_size_ - 1; i >= 0; --i) {
      if (depths_[lhs] - (1 << i) >= depths_[rhs]) {
        lhs = binary_ancestors_[lhs][i];
      }
    }

    if (lhs == rhs) {
      return lhs;
    }

    for (int i = log_tree_size_ - 1; i >= 0; --i) {
      if (binary_ancestors_[lhs][i] != binary_ancestors_[rhs][i]) {
        lhs = binary_ancestors_[lhs][i];
        rhs = binary_ancestors_[rhs][i];
      }
    }

    return binary_ancestors_[lhs][0];
  }

  long long GetDistanceForLCA(long long lca) const {
    return distances_in_bridge_tree_[lca];
  }

 private:
  void DfsForLCA(long long vertex, long long parent, long long depth) {
    depths_[vertex] = depth;
    binary_ancestors_[vertex][0] = parent;

    for (long long i = 1; i < log_tree_size_; ++i) {
      long long prev_ancestor = binary_ancestors_[vertex][i - 1];
      binary_ancestors_[vertex][i] = binary_ancestors_[prev_ancestor][i - 1];
    }

    for (long long neighbor : tree_.GetNeighbors(vertex)) {
      if (neighbor != parent) {
        DfsForLCA(neighbor, vertex, depth + 1);
      }
    }
  }

  const Graph& tree_;
  long long root_;
  const Vector1D<long long>& distances_in_bridge_tree_;
  Vector1D<long long> depths_;
  Vector2D<long long> binary_ancestors_;
  long long log_tree_size_;
};

class GameSolver {
 public:
  GameSolver(long long vertex_numbers, long long edge_numbers,
             long long finish_vertex)
      : edge_numbers_(edge_numbers),
        finish_vertex_(finish_vertex),
        original_graph_(vertex_numbers, edge_numbers) {}

  void AddRoad(long long lhs, long long rhs) {
    original_graph_.AddEdge(lhs, rhs);
  }

  long long GetEdgeCount() { return edge_numbers_; }

  void Preprocess() {
    original_graph_.CalculateDistances(finish_vertex_);

    bridge_finder_ = std::make_unique<BridgeFinder>(original_graph_);
    bridge_finder_->FindBridges();

    component_finder_ =
        std::make_unique<ComponentFinder>(original_graph_, *bridge_finder_);
    component_finder_->FindComponents();

    bridge_tree_ =
        std::make_unique<BridgeTree>(original_graph_, *component_finder_);
    bridge_tree_->BuildBridgeTree(finish_vertex_);

    long long root_component =
        component_finder_->GetComponentId(finish_vertex_);
    lca_finder_ = std::make_unique<LcaFinder>(
        bridge_tree_->GetTreeGraph(), root_component,
        bridge_tree_->GetDistancesInBridgeTree());
  }

  long long AnswerQuery(long long lhs, long long rhs) {
    long long component_lhs = component_finder_->GetComponentId(lhs);
    long long component_rhs = component_finder_->GetComponentId(rhs);
    long long lca = lca_finder_->FindLCA(component_lhs, component_rhs);
    return lca_finder_->GetDistanceForLCA(lca);
  }

 private:
  long long edge_numbers_;
  long long finish_vertex_;

  Graph original_graph_;
  std::unique_ptr<BridgeFinder> bridge_finder_;
  std::unique_ptr<ComponentFinder> component_finder_;
  std::unique_ptr<BridgeTree> bridge_tree_;
  std::unique_ptr<LcaFinder> lca_finder_;
};


std::istream& operator>>(std::istream& stream, GameSolver& solver) {
  for (long long i = 0; i < solver.GetEdgeCount(); ++i) {
    long long lhs;
    long long rhs;
    stream >> lhs >> rhs;
    solver.AddRoad(lhs, rhs);
  }

  return stream;
}


int main() {
  long long vertex_numbers = 0;
  long long edge_numbers = 0;
  long long finish_vertex = 0;
  std::cin >> vertex_numbers >> edge_numbers >> finish_vertex;

  GameSolver solver(vertex_numbers, edge_numbers, finish_vertex);
  std::cin >> solver;

  solver.Preprocess();

  long long query_numbers = 0;
  std::cin >> query_numbers;
  for (long long i = 0; i < query_numbers; ++i) {
    long long lhs;
    long long rhs;
    std::cin >> lhs >> rhs;
    std::cout << solver.AnswerQuery(lhs, rhs) << "\n";
  }

  return 0;
}

