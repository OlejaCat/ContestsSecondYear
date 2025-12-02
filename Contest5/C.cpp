#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>
#include <algorithm>

template <typename T>
using Vector1D = std::vector<T>;

template <typename T>
using Vector2D = Vector1D<Vector1D<T>>;

const long long kInfDistance = std::numeric_limits<long long>::max();

class InterestingGraph {
 private:
  const long long kZeroDistance = 0;

 public:
  struct EdgeInfo {
    long long end;
    long long cost;
    long long time;
  };

  using NodeType = EdgeInfo;

  InterestingGraph() = default;

  InterestingGraph(long long vertex_numbers, long long edge_numbers)
      : vertex_numbers_(vertex_numbers),
        edge_numbers_(edge_numbers),
        graph_(vertex_numbers + 1) {}

  InterestingGraph(const InterestingGraph& other)
      : vertex_numbers_(other.vertex_numbers_),
        edge_numbers_(other.edge_numbers_),
        graph_(other.graph_) {}

  InterestingGraph& operator=(const InterestingGraph& other) {
    if (this != &other) {
      vertex_numbers_ = other.vertex_numbers_;
      edge_numbers_ = other.edge_numbers_;
      graph_ = other.graph_;
    }
    return *this;
  }

  void AddEdge(long long start, long long end, long long cost, long long time) {
    assert(start <= vertex_numbers_);
    assert(end <= vertex_numbers_);

    graph_[start].push_back({end, cost, time});
    graph_[end].push_back({start, cost, time});
  }

  Vector1D<EdgeInfo>& GetNeighbors(long long vertex) {
    assert(vertex <= vertex_numbers_);
    return graph_[vertex];
  }

  const Vector1D<EdgeInfo>& GetNeighbors(long long vertex) const {
    assert(vertex <= vertex_numbers_);
    return graph_[vertex];
  }

  long long GetVertexNumber() const { return vertex_numbers_; }

  long long GetEdgesNumber() const { return edge_numbers_; }

 private:
  long long vertex_numbers_;
  long long edge_numbers_;
  Vector2D<EdgeInfo> graph_;
};

class TimeDijkstra {
 public:
  struct State {
    long long vertex;
    long long time;
    long long cost;
    
    bool operator>(const State& other) const {
      return cost > other.cost;
    }
  };

  using PriorityQueue = std::priority_queue<State, Vector1D<State>, std::greater<State>>;

  TimeDijkstra(const InterestingGraph& graph, long long max_time)
      : graph_(graph), max_time_(max_time) {}

  void Calculate(long long start_vertex, long long target_vertex) {
    start_vertex_ = start_vertex;
    target_vertex_ = target_vertex;
    
    Initialize();
    RunDijkstra();
  }

  bool IsTargetReachable() const {
    long long min_cost = kInfDistance;
    for (long long time = 0; time <= max_time_; ++time) {
      if (cost_[target_vertex_][time] < min_cost) {
        min_cost = cost_[target_vertex_][time];
      }
    }
    return min_cost != kInfDistance;
  }

  long long GetMinCost() const {
    long long min_cost = kInfDistance;
    for (long long time = 0; time <= max_time_; ++time) {
      if (cost_[target_vertex_][time] < min_cost) {
        min_cost = cost_[target_vertex_][time];
      }
    }
    return min_cost;
  }

  Vector1D<long long> RestorePath() const {
    Vector1D<long long> path;
    
    if (!IsTargetReachable()) {
      return path;
    }

    long long optimal_time = -1;
    long long min_cost = kInfDistance;
    for (long long time = 0; time <= max_time_; ++time) {
      if (cost_[target_vertex_][time] < min_cost) {
        min_cost = cost_[target_vertex_][time];
        optimal_time = time;
      }
    }

    long long current_vertex = target_vertex_;
    long long current_time = optimal_time;
    
    while (current_vertex != -1) {
      path.push_back(current_vertex);
      long long prev_vertex = previous_vertex_[current_vertex][current_time];
      long long prev_time = previous_time_[current_vertex][current_time];
      current_vertex = prev_vertex;
      current_time = prev_time;
    }

    std::reverse(path.begin(), path.end());
    return path;
  }

 private:
  void Initialize() {
    long long vertex_count = graph_.GetVertexNumber();
    cost_.assign(vertex_count + 1, Vector1D<long long>(max_time_ + 1, kInfDistance));
    previous_vertex_.assign(vertex_count + 1, Vector1D<long long>(max_time_ + 1, -1));
    previous_time_.assign(vertex_count + 1, Vector1D<long long>(max_time_ + 1, -1));
    
    cost_[start_vertex_][0] = 0;
  }

  void RunDijkstra() {
    PriorityQueue pq;
    pq.push({start_vertex_, 0, 0});

    while (!pq.empty()) {
      State current = pq.top();
      pq.pop();

      long long current_vertex = current.vertex;
      long long current_time = current.time;
      long long current_cost = current.cost;

      if (current_cost != cost_[current_vertex][current_time]) {
        continue;
      }

      for (const auto& edge : graph_.GetNeighbors(current_vertex)) {
        long long neighbor = edge.end;
        long long edge_cost = edge.cost;
        long long edge_time = edge.time;
        
        long long new_time = current_time + edge_time;
        long long new_cost = current_cost + edge_cost;

        if (new_time <= max_time_ && new_cost < cost_[neighbor][new_time]) {
          cost_[neighbor][new_time] = new_cost;
          previous_vertex_[neighbor][new_time] = current_vertex;
          previous_time_[neighbor][new_time] = current_time;
          pq.push({neighbor, new_time, new_cost});
        }
      }
    }
  }

 private:
  InterestingGraph graph_;
  long long max_time_;
  long long start_vertex_;
  long long target_vertex_;
  Vector2D<long long> cost_;
  Vector2D<long long> previous_vertex_;
  Vector2D<long long> previous_time_;
};

void WeCanFindHim();
void InitializeGraph(InterestingGraph& graph);

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  WeCanFindHim();

  return 0;
}

void WeCanFindHim() {
  long long vertex_numbers;
  long long edge_numbers;
  long long max_time;
  std::cin >> vertex_numbers >> edge_numbers >> max_time;

  InterestingGraph graph(vertex_numbers, edge_numbers);
  InitializeGraph(graph);

  TimeDijkstra amogus_finder(graph, max_time);
  amogus_finder.Calculate(1, vertex_numbers);

  if (!amogus_finder.IsTargetReachable()) {
    std::cout << -1 << "\n";
    return;
  }

  long long min_cost = amogus_finder.GetMinCost();
  Vector1D<long long> path = amogus_finder.RestorePath();

  std::cout << min_cost << "\n";
  std::cout << path.size() << "\n";
  for (size_t i = 0; i < path.size(); i++) {
    if (i > 0)
      std::cout << " ";
    std::cout << path[i];
  }
  std::cout << "\n";
}

void InitializeGraph(InterestingGraph& graph) {
  for (long long i = 0; i < graph.GetEdgesNumber(); i++) {
    long long lhs;
    long long rhs;
    long long cost;
    long long time;
    std::cin >> lhs >> rhs >> cost >> time;
    graph.AddEdge(lhs, rhs, cost, time);
  }
}
