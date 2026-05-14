#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <map>
#include <vector>

class GeneralizedSuffixTree {
 private:
  using StateId = std::size_t;
  using Index = std::size_t;
  using StringId = std::size_t;

  static const StateId kRootState = 0;
  static const Index kOpenInterval = std::numeric_limits<Index>::max();
  static const std::size_t kNoParent = std::numeric_limits<std::size_t>::max();

  struct StringCut {
    StringId string_id;
    Index start;
    Index end;

    std::size_t Length() const { return IsEmpty() ? 0 : end - start + 1; }

    bool IsEmpty() const { return end == kOpenInterval || end < start; }
  };

  static constexpr StringCut kEmptyCut = {0, 0, 0};

  struct Edge {
    StringCut cut;
    StateId target_state;
  };

  struct StateType {
    std::map<char, Edge> transitions;
    StateId suffix_link = kRootState;
  };

  struct ActivePoint {
    StateId state;
    Index path_start;
  };

  struct SplitResult {
    bool is_endpoint;
    StateId state_to_extend;
  };

 public:
  GeneralizedSuffixTree() { CreateState(); }

  void Build(std::string_view text, StringId string_id) {
    if (texts_.size() <= string_id) {
      texts_.resize(string_id + 1);
    }

    texts_[string_id] = std::string(text);
    current_string_id_ = string_id;

    ActivePoint active_point = {kRootState, 0};

    for (Index i = 0; i < text.length(); ++i) {
      active_point = Update(active_point.state, active_point.path_start, i);
      active_point = GetCanonicalPoint(active_point.state,
                                       {string_id, active_point.path_start, i}, i);
    }

    for (auto& state : states_) {
      for (auto& [_, edge] : state.transitions) {
        if (edge.cut.string_id == string_id && edge.cut.end == kOpenInterval) {
          edge.cut.end = text.length() - 1;
        }
      }
    }
  }

  void PrintTree() const {
    std::cout << states_.size() << "\n";

    std::size_t id_counter = 0;
    DfsOutput(kRootState, kNoParent, kEmptyCut, id_counter);
  }

 private:
  void DfsOutput(StateId current_state_id, std::size_t parent_logical_id,
                 StringCut incoming_cut, std::size_t& next_id_generator) const {
    std::size_t state_id = next_id_generator;
    ++next_id_generator;

    if (parent_logical_id != kNoParent) {
      Index actual_end_index = 0;
      if (incoming_cut.end == kOpenInterval) {
        actual_end_index = texts_[incoming_cut.string_id].length();
      } else {
        actual_end_index = incoming_cut.end + 1;
      }

      std::cout << parent_logical_id << " " << incoming_cut.string_id << " "
                << incoming_cut.start << " " << actual_end_index << "\n";
    }

    for (const auto& [_, edge] : states_[current_state_id].transitions) {
      DfsOutput(edge.target_state, state_id, edge.cut, next_id_generator);
    }
  }

  ActivePoint Update(StateId active_state, Index path_start, Index current_index) {
    StateId previous_split_state = kRootState;
    char next_symbol = texts_[current_string_id_][current_index];

    StringCut path = {current_string_id_, path_start, current_index - 1};
    SplitResult split = CheckTransitionAndSplit(active_state, path, next_symbol, current_index);

    while (!split.is_endpoint) {
      SetTransition(split.state_to_extend,
                    {current_string_id_, current_index, kOpenInterval},
                    CreateState());

      LinkSuffix(previous_split_state, split.state_to_extend);
      GoSuffixLink(active_state, path_start);

      path = {current_string_id_, path_start, current_index - 1};
      
      ActivePoint canonized_point = GetCanonicalPoint(active_state, path, current_index);
      active_state = canonized_point.state;
      path_start = canonized_point.path_start;

      path = {current_string_id_, path_start, current_index - 1};
      split = CheckTransitionAndSplit(active_state, path, next_symbol, current_index);
    }

    LinkSuffix(previous_split_state, active_state);
    return {active_state, path_start};
  }

  void GoSuffixLink(StateId& active_state, Index& path_start) {
    if (active_state == kRootState) {
      ++path_start;
    }
    active_state = states_[active_state].suffix_link;
  }

  void LinkSuffix(StateId& previous, StateId& current) {
    if (previous != kRootState) {
      states_[previous].suffix_link = current;
    }
    previous = current;
  }

  void SetTransition(StateId source, StringCut cut, StateId target) {
    char edge_start_symbol = texts_[cut.string_id][cut.start];
    states_[source].transitions[edge_start_symbol] = {cut, target};
  }

  ActivePoint GetCanonicalPoint(StateId active_state, StringCut path, Index current_index) {
    if (path.IsEmpty()) {
      return {active_state, path.start};
    }

    auto iterator = states_[active_state].transitions.find(texts_[path.string_id][path.start]);
    if (iterator == states_[active_state].transitions.end()) {
      return {active_state, path.start};
    }

    Edge current_edge = iterator->second;
    Index edge_length = GetRealLength(current_edge.cut, current_index);

    while (edge_length <= path.Length()) {
      path.start += edge_length;
      active_state = current_edge.target_state;
      if (path.IsEmpty()) break;

      iterator = states_[active_state].transitions.find(texts_[path.string_id][path.start]);
      if (iterator == states_[active_state].transitions.end()) break;

      current_edge = iterator->second;
      edge_length = GetRealLength(current_edge.cut, current_index);
    }
    return {active_state, path.start};
  }

  SplitResult CheckTransitionAndSplit(StateId state, StringCut path, char next_symbol, Index current_index) {
    if (path.IsEmpty()) {
      return {states_[state].transitions.count(next_symbol) > 0, state};
    }

    auto iterator = states_[state].transitions.find(texts_[path.string_id][path.start]);
    Edge edge = iterator->second;
    Index split_point = edge.cut.start + path.Length();

    if (next_symbol == texts_[edge.cut.string_id][split_point]) {
      return {true, state};
    }

    StateId new_state = CreateState();
    SetTransition(new_state, {edge.cut.string_id, split_point, edge.cut.end}, edge.target_state);
    SetTransition(state, {edge.cut.string_id, edge.cut.start, split_point - 1}, new_state);

    return {false, new_state};
  }

  StateId CreateState() {
    states_.emplace_back();
    return states_.size() - 1;
  }

  std::size_t GetRealLength(StringCut cut, Index current_index) const {
    if (cut.end == kOpenInterval) {
      return current_index - cut.start + 1;
    }
    return cut.end - cut.start + 1;
  }

  std::vector<std::string> texts_;
  StringId current_string_id_ = 0;
  std::vector<StateType> states_;
};

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::string first;
  std::string second;
  std::cin >> first >> second;

  GeneralizedSuffixTree tree;
  tree.Build(first, 0);
  tree.Build(second, 1);

  tree.PrintTree();
  return 0;
}
