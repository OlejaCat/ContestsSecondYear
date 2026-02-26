#include <algorithm>
#include <array>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

using NodeReference = int;

struct Mapping {
  static const long long kAlphabetPower = 3;
  static NodeReference CharToIndex(char symbol) {
    return symbol - kStartSymbol;
  }

 private:
  static const char kStartSymbol = 'a';
};

class Node {
 public:
  Node() { children_.fill(-1); }

  NodeReference& GetNext(long long index) { return children_[index]; }
  const NodeReference& GetNext(long long index) const {
    return children_[index];
  }

  void AddWord() { ++word_count_; }
  long long WordCount() const { return word_count_; }

  void SetTotalWordCount(long long count) { total_word_count_ = count; }
  long long TotalWordCount() const { return total_word_count_; }

  NodeReference& Link() { return link_; }

 private:
  std::array<NodeReference, Mapping::kAlphabetPower> children_;
  NodeReference link_ = 0;
  long long word_count_ = 0;
  long long total_word_count_ = 0;
};

template <typename Mapping, typename NodeType>
class Trie {
 public:
  Trie() { nodes_.emplace_back(); }

  void Insert(const std::string& str) {
    NodeReference current_node = 0;
    for (char symbol : str) {
      long long index = Mapping::CharToIndex(symbol);
      if (nodes_[current_node].GetNext(index) == -1) {
        nodes_[current_node].GetNext(index) = static_cast<int>(nodes_.size());
        nodes_.emplace_back();
      }
      current_node = nodes_[current_node].GetNext(index);
    }
    nodes_[current_node].AddWord();
  }

  void Reset(size_t expected_nodes = 0) {
    nodes_.clear();
    if (expected_nodes > 0) {
      nodes_.reserve(expected_nodes);
    }
    nodes_.emplace_back();
  }

 protected:
  NodeType& GetNode(NodeReference node_reference) {
    return nodes_[node_reference];
  }

  const NodeType& GetNode(NodeReference node_reference) const {
    return nodes_[node_reference];
  }

  std::vector<NodeType> nodes_;
};

template <typename Mapping, typename NodeType>
class AhoCorasickAutomate : public Trie<Mapping, NodeType> {
 public:
  using Trie<Mapping, NodeType>::GetNode;

  void Build() {
    std::queue<NodeReference> examined_nodes;
    InitializeRoot(examined_nodes);

    while (!examined_nodes.empty()) {
      NodeReference current_node = examined_nodes.front();
      examined_nodes.pop();

      SetTotalWordCountOnNode(current_node);

      for (long long i = 0; i < Mapping::kAlphabetPower; ++i) {
        NodeReference& target = GetNode(current_node).GetNext(i);
        NodeReference link_target =
            GetNode(GetNode(current_node).Link()).GetNext(i);

        if (target != -1) {
          GetNode(target).Link() = link_target;
          examined_nodes.push(target);
        } else {
          target = link_target;
        }
      }
    }
  }

  long long GetTotalOccurancies(const std::string& text) const {
    long long result = 0;
    NodeReference current_node = 0;
    for (char symbol : text) {
      current_node =
          GetNode(current_node).GetNext(Mapping::CharToIndex(symbol));
      result += GetNode(current_node).TotalWordCount();
    }
    return result;
  }

 private:
  void InitializeRoot(std::queue<NodeReference>& examined_nodes) {
    for (long long i = 0; i < Mapping::kAlphabetPower; ++i) {
      NodeReference& next = GetNode(0).GetNext(i);
      if (next != -1) {
        GetNode(next).Link() = 0;
        examined_nodes.push(next);
      } else {
        next = 0;
      }
    }
  }

  void SetTotalWordCountOnNode(NodeReference node) {
    long long node_word_count = GetNode(node).WordCount();
    long long suffix_word_count =
        GetNode(GetNode(node).Link()).TotalWordCount();
    GetNode(node).SetTotalWordCount(node_word_count + suffix_word_count);
  }
};

class DynamicAhoCorasick {
 private:
  struct Block {
    std::vector<std::string> strings;
    AhoCorasickAutomate<Mapping, Node> automate;
  };

 public:
  void Insert(const std::string& str) {
    std::vector<std::string> current_strings = {str};

    while (!blocks_.empty() &&
           blocks_.back().strings.size() == current_strings.size()) {
      for (auto& str : blocks_.back().strings) {
        current_strings.push_back(std::move(str));
      }
      blocks_.pop_back();
    }

    blocks_.emplace_back();
    blocks_.back().strings = current_strings;

    size_t expected_nodes = 1;
    for (const auto& str : current_strings) {
      expected_nodes += str.length();
    }

    blocks_.back().automate.Reset(expected_nodes);
    for (const auto& str : blocks_.back().strings) {
      blocks_.back().automate.Insert(str);
    }
    blocks_.back().automate.Build();
  }

  long long GetTotalOccurancies(const std::string& text) const {
    long long result = 0;
    for (const auto& block : blocks_) {
      result += block.automate.GetTotalOccurancies(text);
    }
    return result;
  }

 private:
  std::vector<Block> blocks_;
};

class RequestsHandler {
   private:
  static const char kAddition = '+';
  static const char kDeletion = '-';
  static const char kOccurancies = '?';
  const std::string kUnknownRequest = "Unknown request type\n";

 public:
  RequestsHandler() = default;

  RequestsHandler(long long requests_count) : requests_count_(requests_count) {}

  void Make(char request_type, std::string& input_string) {
    Shiftstring(input_string);

    switch (request_type) {
      case kAddition:
        HandleAddition(input_string);
        break;

      case kDeletion:
        HandleDeletion(input_string);
        break;

      case kOccurancies:
        HandleOccurancies(input_string);
        break;

      default:
        std::cout << kUnknownRequest;
    }
  }

 private:
  void Shiftstring(std::string& string) const {
    if (string.empty()) {
      return;
    }
    long long shift = last_answer_ % string.length();
    std::rotate(string.begin(), string.begin() + shift, string.end());
  }

  void HandleAddition(const std::string& input_string) {
    if (active_strings_count_[input_string] == 0) {
      active_strings_count_[input_string] = 1;
      positive_forest_.Insert(input_string);
    }
  }

  void HandleDeletion(const std::string& input_string) {
    if (active_strings_count_[input_string] == 1) {
      active_strings_count_[input_string] = 0;
      negative_forest_.Insert(input_string);
    }
  }

  void HandleOccurancies(const std::string& input_string) {
    long long positive_matches =
        positive_forest_.GetTotalOccurancies(input_string);
    long long negative_matches =
        negative_forest_.GetTotalOccurancies(input_string);

    last_answer_ = positive_matches - negative_matches;
    std::cout << last_answer_ << "\n";
  }

  long long requests_count_ = 0;
  long long last_answer_ = 0;

  DynamicAhoCorasick positive_forest_;
  DynamicAhoCorasick negative_forest_;
  std::unordered_map<std::string, int> active_strings_count_;
};

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  long long requests_count = 0;
  std::cin >> requests_count;

  RequestsHandler requests_handler(requests_count);
  for (long long i = 0; i < requests_count; ++i) {
    char request_type;
    std::string input_string;
    std::cin >> request_type >> input_string;

    requests_handler.Make(request_type, input_string);
  }
}
