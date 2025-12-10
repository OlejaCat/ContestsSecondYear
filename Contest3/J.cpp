#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <tuple>
#include <vector>

namespace TagSolver {

enum class Color { WHITE, GRAY, BLACK };

const int kGridSize = 3;

class NodeData {
 public:
  NodeData(Color color = Color::WHITE,
           const std::string& parent_grid_string = "", char move_type = '\0',
           int depth = 0)
      : color_(color),
        parent_grid_string_(parent_grid_string),
        move_type_(move_type),
        depth_(depth) {}

  const Color& GetColor() const { return color_; }
  const std::string& GetParent() const { return parent_grid_string_; }
  const char& GetMove() const { return move_type_; }
  const int& GetDepth() const { return depth_; }

  void SetColor(Color color) { color_ = color; }

 private:
  Color color_;
  std::string parent_grid_string_;
  char move_type_;
  int depth_;
};

class PuzzleSolver {
 private:
  const std::string kSolvedGrid = "123456780";

  static const char kUpMove = 'U';
  static const char kDownMove = 'D';
  static const char kLeftMove = 'L';
  static const char kRightMove = 'R';

  const std::vector<std::tuple<int, int, char>> kAvailibleMoves{
      {-1, +0, kUpMove},
      {+1, +0, kDownMove},
      {+0, -1, kLeftMove},
      {+0, +1, kRightMove}};

 public:
  PuzzleSolver(const std::string& initial_grid_string)
      : initial_grid_(initial_grid_string) {
    if (initial_grid_string != kSolvedGrid) {
      queue_forward_.push(initial_grid_string);
      forward_search_[initial_grid_string] = NodeData(Color::GRAY, "", '\0', 0);

      queue_reverse_.push(kSolvedGrid);
      reverse_search_[kSolvedGrid] = NodeData(Color::GRAY, "", '\0', 0);
    }
  }

  void Solve() {
    if (initial_grid_ == kSolvedGrid) {
      std::cout << 0 << "\n";
      return;
    }

    std::string meeting_key;

    while (!queue_forward_.empty() && !queue_reverse_.empty()) {
      bool found = false;

      if (queue_forward_.size() <= queue_reverse_.size()) {
        found = SearchStep(queue_forward_, forward_search_, reverse_search_,
                           meeting_key);
      } else {
        found = SearchStep(queue_reverse_, reverse_search_, forward_search_,
                           meeting_key);
      }

      if (found) {
        RestorePath(meeting_key);
        return;
      }
    }

    std::cout << -1 << "\n";
  }

 private:
  static int GetEmptyCell(const std::string& grid_string) {
    return static_cast<int>(grid_string.find('0'));
  }

  std::vector<std::pair<std::string, char>> GetTransitions(
      const std::string& current_grid_string) {
    std::vector<std::pair<std::string, char>> transitions;
    int empty_cell_index = GetEmptyCell(current_grid_string);

    int row = empty_cell_index / kGridSize;
    int column = empty_cell_index % kGridSize;

    for (const auto& move : kAvailibleMoves) {
      int new_row = row + std::get<0>(move);
      int new_column = column + std::get<1>(move);
      char move_type = std::get<2>(move);

      if (CheckBorders(new_row) && CheckBorders(new_column)) {
        int new_index = new_row * kGridSize + new_column;

        std::string new_grid_string = current_grid_string;
        std::swap(new_grid_string[empty_cell_index],
                  new_grid_string[new_index]);

        transitions.emplace_back(new_grid_string, move_type);
      }
    }
    return transitions;
  }

  static bool CheckBorders(int index) {
    return index >= 0 && index < kGridSize;
  }

  bool SearchStep(std::queue<std::string>& current_queue,
                  std::map<std::string, NodeData>& current_search,
                  std::map<std::string, NodeData>& opposite_search,
                  std::string& meeting_state) {
    if (current_queue.empty()) {
      return false;
    }

    std::string current_grid = current_queue.front();
    current_queue.pop();

    auto transitions = GetTransitions(current_grid);

    for (const auto& transition : transitions) {
      const std::string& next_grid = transition.first;
      char move_type = transition.second;

      if (current_search.find(next_grid) == current_search.end()) {
        int new_depth = current_search[current_grid].GetDepth() + 1;
        current_search[next_grid] =
            NodeData(Color::GRAY, current_grid, move_type, new_depth);
        current_queue.push(next_grid);
      }

      if (opposite_search.find(next_grid) != opposite_search.end()) {
        meeting_state = next_grid;
        return true;
      }
    }

    current_search[current_grid].SetColor(Color::BLACK);
    return false;
  }

  static char GetTileMove(char zero_move) {
    switch (zero_move) {
      case kDownMove:
        return kUpMove;
      case kUpMove:
        return kDownMove;
      case kLeftMove:
        return kRightMove;
      case kRightMove:
        return kLeftMove;
      default:
        return '\0';
    }
  }

  std::string RestoreForwardPath(const std::string& meeting_state) {
    std::string forward_path;
    std::string current_key = meeting_state;

    while (!forward_search_[current_key].GetParent().empty()) {
      char zero_move = forward_search_[current_key].GetMove();
      forward_path += zero_move;
      current_key = forward_search_[current_key].GetParent();
    }

    std::reverse(forward_path.begin(), forward_path.end());
    return forward_path;
  }

  std::string RestoreReversePath(const std::string& meeting_state) {
    std::string reverse_path;
    std::string current_key = meeting_state;

    while (!reverse_search_[current_key].GetParent().empty()) {
      char zero_move = reverse_search_[current_key].GetMove();
      reverse_path += GetTileMove(zero_move);
      current_key = reverse_search_[current_key].GetParent();
    }

    return reverse_path;
  }

  void RestorePath(const std::string& meeting_state) {
    std::string forward_path = RestoreForwardPath(meeting_state);
    std::string reverse_path = RestoreReversePath(meeting_state);

    int total_depth = forward_search_[meeting_state].GetDepth()
                      + reverse_search_[meeting_state].GetDepth();

    std::cout << total_depth << "\n";
    std::cout << forward_path << reverse_path << "\n";
  }

  std::string initial_grid_;

  std::map<std::string, NodeData> forward_search_;
  std::map<std::string, NodeData> reverse_search_;

  std::queue<std::string> queue_forward_;
  std::queue<std::string> queue_reverse_;
};

}  // namespace TagSolver

int main() {
  std::string initial_grid;
  int value;
  for (int i = 0; i < TagSolver::kGridSize * TagSolver::kGridSize; ++i) {
    std::cin >> value;
    initial_grid += std::to_string(value);
  }

  TagSolver::PuzzleSolver solution(initial_grid);
  solution.Solve();

  return 0;
}
