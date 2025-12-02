#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

template <typename T>
using Vector1D = std::vector<T>;

template <typename T>
using Vector2D = Vector1D<Vector1D<T>>;

const long long kInfDistance = std::numeric_limits<long long>::max();

class WeightedUndirectedGraph {
 private:
  const long long kZeroDistance = 0;

 public:
  using NodeType = std::pair<long long, long long>;

  WeightedUndirectedGraph() = default;

  WeightedUndirectedGraph(long long vertex_numbers, long long edge_numbers)
      : vertex_numbers_(vertex_numbers),
        edge_numbers_(edge_numbers),
        graph_(vertex_numbers + 1) {}

  WeightedUndirectedGraph(const WeightedUndirectedGraph& other)
      : vertex_numbers_(other.vertex_numbers_),
        edge_numbers_(other.edge_numbers_),
        graph_(other.graph_) {}

  WeightedUndirectedGraph& operator=(const WeightedUndirectedGraph& other) {
    if (this != &other) {
      vertex_numbers_ = other.vertex_numbers_;
      edge_numbers_ = other.edge_numbers_;
      graph_ = other.graph_;
    }
    return *this;
  }

  void AddEdge(long long start, long long end, long long weight) {
    assert(start <= vertex_numbers_);
    assert(end <= vertex_numbers_);

    graph_[start].push_back({end, weight});
    graph_[end].push_back({start, weight});
  }

  Vector1D<NodeType>& GetNeighbors(long long vertex) {
    assert(vertex <= vertex_numbers_);

    return graph_[vertex];
  }

  const Vector1D<NodeType>& GetNeighbors(long long vertex) const {
    assert(vertex <= vertex_numbers_);

    return graph_[vertex];
  }

  long long GetVertexNumber() const { return vertex_numbers_; }

  long long GetEdgesNumber() const { return edge_numbers_; }

 private:
  long long vertex_numbers_;
  long long edge_numbers_;
  Vector2D<NodeType> graph_;
};

class DijkstraBase {
 public:
  struct QueueNode {
    long long distance;
    long long vertex;

    bool operator>(const QueueNode& other) const {
      return distance > other.distance;
    }
  };

  using PriorityQueue = std::priority_queue<QueueNode, Vector1D<QueueNode>,
                                            std::greater<QueueNode>>;

  DijkstraBase(const WeightedUndirectedGraph& graph) : graph_(graph) {}

  virtual void Calculate() = 0;

  Vector1D<long long>& GetDistances() { return distances_; }

  const Vector1D<long long>& GetDistances() const { return distances_; }

  long long GetDistanceTo(long long vertex) const {
    assert(vertex <= graph_.GetVertexNumber());

    return distances_[vertex];
  }

 protected:
  WeightedUndirectedGraph graph_;
  Vector1D<long long> distances_;

  virtual bool ShouldUpdateWeight(long long current_vertex, long long neighbor,
                                  long long current_distance,
                                  long long edge_weight) {
    long long new_distance = current_distance + edge_weight;
    return new_distance < distances_[neighbor];
  }

  virtual void Initialize() {
    distances_.assign(graph_.GetVertexNumber() + 1, kInfDistance);
  }

  void RunDijkstra() {
    PriorityQueue vertex_to_process;

    InitializeQueue(vertex_to_process);

    while (!vertex_to_process.empty()) {
      QueueNode node = vertex_to_process.top();
      vertex_to_process.pop();

      long long current_distance = node.distance;
      long long current_vertex = node.vertex;

      if (current_distance != distances_[current_vertex]) {
        continue;
      }

      for (auto [neighbor, edge_weight] : graph_.GetNeighbors(current_vertex)) {
        if (ShouldUpdateWeight(current_vertex, neighbor, current_distance,
                               edge_weight)) {
          distances_[neighbor] = current_distance + edge_weight;
          vertex_to_process.push({distances_[neighbor], neighbor});
        }
      }
    }
  }

  virtual void InitializeQueue(PriorityQueue& vertex_to_process) = 0;
};

class GraphInfection : public DijkstraBase {
 public:
  GraphInfection(const WeightedUndirectedGraph& graph,
                 const Vector1D<long long>& sources)
      : DijkstraBase(graph), sources_(sources) {}

  void Calculate() override {
    Initialize();
    RunDijkstra();
  }

 protected:
  void InitializeQueue(PriorityQueue& vertex_to_process) override {
    for (const auto& source : sources_) {
      distances_[source] = 0;
      vertex_to_process.push({0, source});
    }
  }

 private:
  Vector1D<long long> sources_;
};

class PlayerDijkstra : public DijkstraBase {
 public:
  PlayerDijkstra(const WeightedUndirectedGraph& graph, long long start,
                 const Vector1D<long long>& virus_distances)
      : DijkstraBase(graph), start_(start), virus_distances_(virus_distances) {}

  void Calculate() override {
    Initialize();
    RunDijkstra();
  }

 protected:
  bool ShouldUpdateWeight(long long current_vertex, long long neighbor,
                          long long current_distance,
                          long long edge_weight) override {
    long long new_distance = current_distance + edge_weight;
    return new_distance < virus_distances_[neighbor]
           && new_distance < distances_[neighbor];
  }

  void InitializeQueue(PriorityQueue& vertex_to_process) override {
    distances_[start_] = 0;
    vertex_to_process.push({0, start_});
  }

 private:
  long long start_;
  Vector1D<long long> virus_distances_;
};

void CanYouWin();
void InitializeGraph(WeightedUndirectedGraph& graph);

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(0);

  CanYouWin();

  return 0;
}

void CanYouWin() {
  long long vertex_number;
  long long edge_number;
  long long virus_numbers;
  std::cin >> vertex_number >> edge_number >> virus_numbers;

  Vector1D<long long> infected(virus_numbers);
  for (long long i = 0; i < virus_numbers; ++i) {
    std::cin >> infected[i];
  }

  WeightedUndirectedGraph graph(vertex_number, edge_number);
  InitializeGraph(graph);

  long long start_point, cure_point;
  std::cin >> start_point >> cure_point;

  GraphInfection virus_dijkstra(graph, infected);
  virus_dijkstra.Calculate();
  const auto& virus_distances = virus_dijkstra.GetDistances();

  PlayerDijkstra player_dijkstra(graph, start_point, virus_distances);
  player_dijkstra.Calculate();

  long long result = player_dijkstra.GetDistanceTo(cure_point);
  if (result == kInfDistance) {
    std::cout << -1 << "\n";
  } else {
    std::cout << result << "\n";
  }
}

void InitializeGraph(WeightedUndirectedGraph& graph) {
  for (long long i = 0; i < graph.GetEdgesNumber(); ++i) {
    long long lhs;
    long long rhs;
    long long weight;
    std::cin >> lhs >> rhs >> weight;
    graph.AddEdge(lhs, rhs, weight);
  }
}
