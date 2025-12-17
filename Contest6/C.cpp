#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

template <typename T>
using Vector1D = std::vector<T>;

using VertexNumberType = long long;
using EdgeIndexType = std::size_t;

const VertexNumberType kSourceVertex = 1;

class FlowGraph {
 private:
  class EdgeInfo {
   public:
    explicit EdgeInfo(VertexNumberType neighbor, long long capacity,
                      EdgeIndexType reverse_edge_index)
        : neighbor_(neighbor),
          flow_(0),
          capacity_(capacity),
          reverse_edge_index_(reverse_edge_index) {}

    auto GetResidualCapacity() const { return capacity_ - flow_; }
    VertexNumberType GetNeighbor() const { return neighbor_; }
    EdgeIndexType GetReverseEdgeIndex() const { return reverse_edge_index_; }

    void IncreaseFlow(long long delta) { flow_ += delta; }
    long long GetFlow() const { return flow_; }
    long long GetCapacity() const { return capacity_; }

   private:
    VertexNumberType neighbor_;
    long long flow_;
    long long capacity_;
    EdgeIndexType reverse_edge_index_;
  };

 public:
  FlowGraph() = default;

  FlowGraph(VertexNumberType vertices_number, VertexNumberType edges_number,
            VertexNumberType source, VertexNumberType sink)
      : adjacement_graph_(vertices_number + 1),
        edges_number_(edges_number),
        source_(source),
        sink_(sink) {}

  void AddEdge(VertexNumberType start, VertexNumberType end,
               long long capacity) {
    EdgeIndexType forward_index = edge_info_graph_.size();
    EdgeIndexType reverse_index = edge_info_graph_.size() + 1;

    edge_info_graph_.emplace_back(end, capacity, reverse_index);
    adjacement_graph_[start].push_back(forward_index);

    edge_info_graph_.emplace_back(start, 0, forward_index);
    adjacement_graph_[end].push_back(reverse_index);
  }

  VertexNumberType GetNeighbor(VertexNumberType vertex,
                               EdgeIndexType neighbor_index) const {
    EdgeIndexType edge_idx = adjacement_graph_[vertex][neighbor_index];
    return edge_info_graph_[edge_idx].GetNeighbor();
  }

  long long GetResidualCapacity(VertexNumberType vertex,
                                EdgeIndexType neighbor_index) const {
    EdgeIndexType edge_idx = adjacement_graph_[vertex][neighbor_index];
    return edge_info_graph_[edge_idx].GetResidualCapacity();
  }

  void IncreaseFlow(VertexNumberType vertex, EdgeIndexType neighbor_index,
                    long long flow_increase) {
    EdgeIndexType forward_edge_idx = adjacement_graph_[vertex][neighbor_index];
    edge_info_graph_[forward_edge_idx].IncreaseFlow(flow_increase);

    EdgeIndexType reverse_edge_idx =
        edge_info_graph_[forward_edge_idx].GetReverseEdgeIndex();
    edge_info_graph_[reverse_edge_idx].IncreaseFlow(-flow_increase);
  }

  std::size_t GetNeighborCount(VertexNumberType vertex) const {
    return adjacement_graph_[vertex].size();
  }

  VertexNumberType GetVerticesNumber() const {
    return static_cast<VertexNumberType>(adjacement_graph_.size() - 1);
  }
  VertexNumberType GetEdgesNumber() const { return edges_number_; }

  VertexNumberType GetSource() const { return source_; }
  VertexNumberType GetSink() const { return sink_; }

 private:
  Vector1D<Vector1D<EdgeIndexType>> adjacement_graph_;
  Vector1D<EdgeInfo> edge_info_graph_;

  VertexNumberType edges_number_;
  VertexNumberType source_;
  VertexNumberType sink_;
};

std::istream& operator>>(std::istream& stream, FlowGraph& flow_graph) {
  for (VertexNumberType i = 1; i <= flow_graph.GetEdgesNumber(); ++i) {
    VertexNumberType start;
    VertexNumberType end;
    long long capacity;
    std::cin >> start >> end >> capacity;
    flow_graph.AddEdge(start, end, capacity);
  }

  return stream;
}

class FordFulkerson {
 public:
  FordFulkerson(FlowGraph& flow_graph) : flow_graph_(flow_graph) {
    used_vertices_in_path_.resize(flow_graph_.GetVerticesNumber() + 1, 0);
  }

  long long GetMaxFlow() {
    long long max_flow = 0;
    long long path_flow;

    while ((path_flow = GetMinIncreasingFlow()) > 0) {
      max_flow += path_flow;
      ++timer_;
    }
    return max_flow;
  }

 private:
  long long GetMinIncreasingFlow() {
    return GetMinIncreasingFlowImpl(flow_graph_.GetSource(),
                                    std::numeric_limits<long long>::max());
  }

  long long GetMinIncreasingFlowImpl(VertexNumberType vertex,
                                     long long min_capacity) {
    if (vertex == flow_graph_.GetSink()) {
      return min_capacity;
    }

    used_vertices_in_path_[vertex] = timer_;

    for (EdgeIndexType i = 0; i < flow_graph_.GetNeighborCount(vertex); ++i) {
      VertexNumberType neighbor = flow_graph_.GetNeighbor(vertex, i);
      long long residual_capacity = flow_graph_.GetResidualCapacity(vertex, i);

      if (residual_capacity == 0 ||
          used_vertices_in_path_[neighbor] == timer_) {
        continue;
      }

      long long flow_increase = GetMinIncreasingFlowImpl(
          neighbor, std::min(min_capacity, residual_capacity));

      if (flow_increase > 0) {
        flow_graph_.IncreaseFlow(vertex, i, flow_increase);
        return flow_increase;
      }
    }

    return 0;
  }

  FlowGraph& flow_graph_;
  Vector1D<std::size_t> used_vertices_in_path_;
  std::size_t timer_ = 1;
};

int main() {
  VertexNumberType vertices_number;
  VertexNumberType edges_number;
  std::cin >> vertices_number >> edges_number;

  FlowGraph flow_graph(vertices_number, edges_number, kSourceVertex,
                       vertices_number);
  std::cin >> flow_graph;

  FordFulkerson ford_fulkerson_algorithm(flow_graph);
  std::cout << ford_fulkerson_algorithm.GetMaxFlow() << "\n";

  return 0;
}
