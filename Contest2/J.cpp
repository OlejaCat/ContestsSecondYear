#include <algorithm>
#include <iostream>
#include <vector>

static constexpr long long kMaxValue = 1e12;
static constexpr long long kMinValue = -1e12;

class ConvexHullTrick {
 private:
  struct Line {
    Line(long long coeff_k, long long coeff_b)
        : coeff_k(coeff_k), coeff_b(coeff_b) {}

    long long GetValue(long long dot) const {
      return (coeff_k * dot) + coeff_b;
    }

    friend long long Cross(const Line& lhs, const Line& rhs) {
      long long numerator = rhs.coeff_b - lhs.coeff_b;
      long long denominator = lhs.coeff_k - rhs.coeff_k;
      if (numerator > 0) {
        return (numerator + denominator - 1) / denominator;
      }

      return numerator / denominator;
    }

    long long coeff_k;
    long long coeff_b;
  };

 public:
  ConvexHullTrick() : dots_{kMinValue} {}

  void AddLine(long long coeff_k, long long coeff_b) {
    Line new_line(coeff_k, coeff_b);
    while (!lines_.empty() && CheckUseless(lines_.back(), new_line)) {
      lines_.pop_back();
      dots_.pop_back();
    }

    if (!lines_.empty()) {
      dots_.push_back(Cross(lines_.back(), new_line));
    }
    lines_.push_back(new_line);
  }

  long long GetMin(long long dot) {
    if (lines_.empty()) {
      return kMaxValue;
    }

    auto iterator = std::upper_bound(dots_.begin(), dots_.end(), dot);
    auto position = iterator - dots_.begin();
    return lines_[position - 1].GetValue(dot);
  }

 private:
  bool CheckUseless(const Line& lhs, const Line& rhs) {
    return lhs.GetValue(dots_.back()) > rhs.GetValue(dots_.back());
  }

  std::vector<Line> lines_;
  std::vector<long long> dots_;
};

long long GetJumpsMin(std::size_t number_of_stones,
                      const std::vector<long long>& stones_heights,
                      long long cost_const);

int main() {
  std::size_t number_of_stones = 0;
  long long cost_const = 0;

  std::cin >> number_of_stones >> cost_const;
  std::vector<long long> stones_heights(number_of_stones + 1, 0);
  for (std::size_t i = 1; i <= number_of_stones; ++i) {
    std::cin >> stones_heights[i];
  }

  std::cout << GetJumpsMin(number_of_stones, stones_heights, cost_const)
            << "\n";
  return 0;
}

long long GetJumpsMin(std::size_t number_of_stones,
                      const std::vector<long long>& stones_heights,
                      long long cost_const) {
  std::vector<long long> min_costs(number_of_stones + 1, kMaxValue);
  min_costs[1] = 0;
  ConvexHullTrick convex_hull_trick;

  convex_hull_trick.AddLine(-2 * stones_heights[1],
                            stones_heights[1] * stones_heights[1]);
  for (std::size_t i = 2; i <= number_of_stones; ++i) {
    long long current_min_value = convex_hull_trick.GetMin(stones_heights[i]);
    min_costs[i] =
        stones_heights[i] * stones_heights[i] + cost_const + current_min_value;
    convex_hull_trick.AddLine(
        -2 * stones_heights[i],
        min_costs[i] + (stones_heights[i] * stones_heights[i]));
  }
  return min_costs[number_of_stones];
}

