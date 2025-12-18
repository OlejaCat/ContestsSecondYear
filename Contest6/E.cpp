#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

template <typename T>
using Vector1D = std::vector<T>;

using VertexNumberType = long long;
using EdgeIndexType = std::size_t;

const long long kReverseEdgeId = -1;

class FlowGraph {
 private:
  class EdgeInfo {
   public:
    explicit EdgeInfo(VertexNumberType neighbor, long long capacity,
                      EdgeIndexType reverse_edge_index, long long original_id)
        : neighbor_(neighbor),
          flow_(0),
          capacity_(capacity),
          reverse_edge_index_(reverse_edge_index),
          original_edge_id_(original_id) {}

    auto GetResidualCapacity() const { return capacity_ - flow_; }
    VertexNumberType GetNeighbor() const { return neighbor_; }
    EdgeIndexType GetReverseEdgeIndex() const { return reverse_edge_index_; }

    void IncreaseFlow(long long delta) { flow_ += delta; }
    long long GetFlow() const { return flow_; }
    long long GetCapacity() const { return capacity_; }
    long long GetOriginalId() const { return original_edge_id_; }

   private:
    VertexNumberType neighbor_;
    long long flow_;
    long long capacity_;
    EdgeIndexType reverse_edge_index_;

    long long original_edge_id_ = kReverseEdgeId;
  };

 public:
  FlowGraph() = default;

  FlowGraph(VertexNumberType vertices_number, VertexNumberType edges_number,
            VertexNumberType source, VertexNumberType sink)
      : adjacency_list_(vertices_number + 1),
        edges_number_(edges_number),
        source_(source),
        sink_(sink) {}

  void AddEdge(VertexNumberType start, VertexNumberType end,
               long long capacity) {
    EdgeIndexType forward_index = edge_info_vector_.size();
    EdgeIndexType reverse_index = edge_info_vector_.size() + 1;
    long long current_id = next_original_edge_id_++;

    edge_info_vector_.emplace_back(end, capacity, reverse_index, current_id);
    adjacency_list_[start].push_back(forward_index);

    edge_info_vector_.emplace_back(start, 0, forward_index, -1);
    adjacency_list_[end].push_back(reverse_index);
  }

  VertexNumberType GetNeighbor(VertexNumberType vertex,
                               EdgeIndexType neighbor_index) const {
    EdgeIndexType edge_idx = adjacency_list_[vertex][neighbor_index];
    return edge_info_vector_[edge_idx].GetNeighbor();
  }

  long long GetResidualCapacity(VertexNumberType vertex,
                                EdgeIndexType neighbor_index) const {
    EdgeIndexType edge_idx = adjacency_list_[vertex][neighbor_index];
    return edge_info_vector_[edge_idx].GetResidualCapacity();
  }

  void IncreaseFlow(VertexNumberType vertex, EdgeIndexType neighbor_index,
                    long long flow_increase) {
    EdgeIndexType forward_edge_idx = adjacency_list_[vertex][neighbor_index];
    edge_info_vector_[forward_edge_idx].IncreaseFlow(flow_increase);

    EdgeIndexType reverse_edge_idx =
        edge_info_vector_[forward_edge_idx].GetReverseEdgeIndex();
    edge_info_vector_[reverse_edge_idx].IncreaseFlow(-flow_increase);
  }

  long long GetActualFlow(VertexNumberType vertex,
                          EdgeIndexType neighbor_index) const {
    EdgeIndexType edge_idx = adjacency_list_[vertex][neighbor_index];
    return edge_info_vector_[edge_idx].GetFlow();
  }

  long long GetCapacity(VertexNumberType vertex,
                        EdgeIndexType neighbor_index) const {
    EdgeIndexType edge_idx = adjacency_list_[vertex][neighbor_index];
    return edge_info_vector_[edge_idx].GetCapacity();
  }

  std::size_t GetNeighborCount(VertexNumberType vertex) const {
    return adjacency_list_[vertex].size();
  }

  VertexNumberType GetVerticesNumber() const {
    return static_cast<VertexNumberType>(adjacency_list_.size() - 1);
  }
  VertexNumberType GetEdgesNumber() const { return edges_number_; }

  VertexNumberType GetSource() const { return source_; }
  VertexNumberType GetSink() const { return sink_; }

  void PrintFlows() const {
    Vector1D<long long> results(next_original_edge_id_, 0);

    for (const auto& edge : edge_info_vector_) {
      long long original_id = edge.GetOriginalId();
      if (original_id > 0) {
        results[original_id] = edge.GetFlow();
      }
    }

    for (long long id = 1; id < next_original_edge_id_; ++id) {
      std::cout << results[id] << "\n";
    }
  }

 private:
  Vector1D<Vector1D<EdgeIndexType>> adjacency_list_;
  Vector1D<EdgeInfo> edge_info_vector_;

  VertexNumberType edges_number_;
  VertexNumberType source_;
  VertexNumberType sink_;

  long long next_original_edge_id_ = 1;
};

std::istream& operator>>(std::istream& stream, FlowGraph& flow_graph) {
  for (VertexNumberType i = 1; i <= flow_graph.GetEdgesNumber(); ++i) {
    VertexNumberType start;
    VertexNumberType end;
    long long capacity;
    stream >> start >> end >> capacity;
    flow_graph.AddEdge(start, end, capacity);
  }

  return stream;
}

class Dinic {
 public:
  Dinic(FlowGraph& flow_graph)
      : flow_graph_(flow_graph),
        level_vector_(flow_graph.GetVerticesNumber() + 1),
        pointer_vector_(flow_graph.GetVerticesNumber() + 1) {}

  long long GetMaxFlow() {
    long long max_flow = 0;

    while (BuildLevelGraph()) {
      std::fill(pointer_vector_.begin(), pointer_vector_.end(), 0);
      long long path_flow;

      while ((path_flow = FindBlockingFlow(
                  flow_graph_.GetSource(),
                  std::numeric_limits<long long>::max())) > 0) {
        max_flow += path_flow;
      }
    }
    return max_flow;
  }

 private:
  bool BuildLevelGraph() {
    std::fill(level_vector_.begin(), level_vector_.end(), -1);

    VertexNumberType source = flow_graph_.GetSource();
    VertexNumberType sink = flow_graph_.GetSink();

    level_vector_[source] = 0;
    std::queue<VertexNumberType> vertex_queue;
    vertex_queue.push(source);

    while (!vertex_queue.empty()) {
      VertexNumberType current_vertex = vertex_queue.front();
      vertex_queue.pop();

      for (EdgeIndexType i = 0;
           i < flow_graph_.GetNeighborCount(current_vertex); ++i) {
        VertexNumberType neighbor = flow_graph_.GetNeighbor(current_vertex, i);
        long long capacity = flow_graph_.GetResidualCapacity(current_vertex, i);

        if (capacity > 0 && level_vector_[neighbor] == -1) {
          level_vector_[neighbor] = level_vector_[current_vertex] + 1;
          vertex_queue.push(neighbor);
        }
      }
    }

    return level_vector_[sink] != -1;
  }

  long long FindBlockingFlow(VertexNumberType vertex, long long flow_limit) {
    if (flow_limit == 0) {
      return 0;
    }

    VertexNumberType sink = flow_graph_.GetSink();
    if (vertex == sink) {
      return flow_limit;
    }

    for (; pointer_vector_[vertex] < flow_graph_.GetNeighborCount(vertex);
         ++pointer_vector_[vertex]) {
      EdgeIndexType idx = pointer_vector_[vertex];
      VertexNumberType neighbor = flow_graph_.GetNeighbor(vertex, idx);
      long long capacity = flow_graph_.GetResidualCapacity(vertex, idx);

      if (level_vector_[neighbor] == level_vector_[vertex] + 1 &&
          capacity > 0) {
        long long pushed_flow =
            FindBlockingFlow(neighbor, std::min(flow_limit, capacity));

        if (pushed_flow > 0) {
          flow_graph_.IncreaseFlow(vertex, idx, pushed_flow);
          return pushed_flow;
        }
      }
    }

    return 0;
  }

  FlowGraph& flow_graph_;
  Vector1D<long long> level_vector_;
  Vector1D<EdgeIndexType> pointer_vector_;
};

struct TranstedMacProblemInput {
  long long number_of_resellers;
  long long number_of_trust_pairs;
  Vector1D<long long> initial_mac_count;
  Vector1D<std::pair<long long, long long>> trust_pairs;
  long long total_initial_mac_count = 0;
  long long maximum_initial_count = 0;
};

struct FlowGraphParameters {
  FlowGraph flow_graph;
  long long requires_flow = 0;
};

std::istream& operator>>(std::istream& stream, TranstedMacProblemInput& input) {
  stream >> input.number_of_resellers >> input.number_of_trust_pairs;
  input.initial_mac_count.resize(input.number_of_resellers);
  for (long long i = 0; i < input.number_of_resellers; ++i) {
    long long mac_count;
    stream >> mac_count;
    input.initial_mac_count[i] = mac_count;

    input.total_initial_mac_count += mac_count;
    input.maximum_initial_count =
        std::max(input.maximum_initial_count, mac_count);
  }

  for (long long i = 0; i < input.number_of_trust_pairs; ++i) {
    long long start;
    long long end;
    stream >> start >> end;
    input.trust_pairs.emplace_back(start, end);
  }

  return stream;
}

FlowGraphParameters InitializeFlowGraph(long long max_load_limit,
                                        const TranstedMacProblemInput& input) {
  VertexNumberType number_of_vertices = input.number_of_resellers + 2;
  VertexNumberType source_vertex = input.number_of_resellers + 1;
  VertexNumberType sink = input.number_of_resellers + 2;
  VertexNumberType edge_number =
      input.number_of_resellers * 2 + input.number_of_trust_pairs;

  FlowGraph flow_graph(number_of_vertices, edge_number, source_vertex, sink);
  long long required_excess_flow = 0;

  for (long long i = 0; i < input.number_of_resellers; ++i) {
    long long current_mac_count = input.initial_mac_count[i];
    VertexNumberType reseller_vertex = i + 1;

    if (current_mac_count > max_load_limit) {
      long long excess_to_transfer = current_mac_count - max_load_limit;
      flow_graph.AddEdge(source_vertex, reseller_vertex, excess_to_transfer);
      required_excess_flow += excess_to_transfer;
    }

    if (current_mac_count < max_load_limit) {
      long long capacity_to_receive = max_load_limit - current_mac_count;
      flow_graph.AddEdge(reseller_vertex, sink, capacity_to_receive);
    }
  }

  long long infinit_capacity = required_excess_flow + 1;
  for (const auto& pair : input.trust_pairs) {
    flow_graph.AddEdge(pair.first, pair.second, infinit_capacity);
  }

  return {std::move(flow_graph), required_excess_flow};
}

bool CheckFlowCondition(FlowGraphParameters& data) {
  if (data.requires_flow == 0) {
    return true;
  }

  Dinic dinic(data.flow_graph);
  long long max_flow_achieved = dinic.GetMaxFlow();

  return max_flow_achieved == data.requires_flow;
}

bool IsPossibleToAchieveLoad(long long max_load_limit,
                             const TranstedMacProblemInput& input) {
  FlowGraphParameters flow_graph_data =
      InitializeFlowGraph(max_load_limit, input);
  return CheckFlowCondition(flow_graph_data);
}

long long FindMinimumMaxLoad(const TranstedMacProblemInput& input) {
  long long total_macs = input.total_initial_mac_count;
  long long num_resellers = input.number_of_resellers;
  long long calculated_min_load =
      (total_macs + num_resellers - 1) / num_resellers;
  long long low = calculated_min_load;
  long long high = input.total_initial_mac_count;

  long long min_max_load = high;

  while (low <= high) {
    long long current_load = low + ((high - low) / 2);

    if (IsPossibleToAchieveLoad(current_load, input)) {
      min_max_load = current_load;
      high = current_load - 1;
    } else {
      low = current_load + 1;
    }
  }

  return min_max_load;
}

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(NULL);

  TranstedMacProblemInput initial_data;
  std::cin >> initial_data;

  std::cout << FindMinimumMaxLoad(initial_data) << "\n";

  return 0;
}
