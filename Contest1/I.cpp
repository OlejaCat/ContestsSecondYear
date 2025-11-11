#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>



void CalculateLongestSequences(const std::vector<int>& gangsters,
                               std::vector<int>& greater_numbers,
                               std::vector<int>& greater_path_choices,
                               std::vector<int>& lesser_numbers,
                               std::vector<int>& lesser_path_choices);

std::tuple<int, int, bool> FindLongestSequenc(
        const std::vector<int>& greater_numbers,
        const std::vector<int>& lesser_numbers);

std::vector<int> RestoreBestSequence(int start_index, bool increasing_order,
                                     std::vector<int>& greater_path_choices,
                                     std::vector<int>& lesser_path_choices,
                                     const std::vector<int>& gangsters);

void Solve();


int main() {
  Solve();

  return 0;
}


void Solve() {
  int number_of_districts = 0;
  std::cin >> number_of_districts;
  std::vector<int> numbers_of_gangsters(number_of_districts);
  for (int i = 0; i < number_of_districts; ++i) {
    std::cin >> numbers_of_gangsters[i];
  }

  std::vector<int> greater_numbers(number_of_districts, 1);
  std::vector<int> greater_path_choices(number_of_districts, -1);
  std::vector<int> lesser_numbers(number_of_districts, 1);
  std::vector<int> lesser_path_choices(number_of_districts, -1);

  CalculateLongestSequences(numbers_of_gangsters, greater_numbers,
                            greater_path_choices, lesser_numbers,
                            lesser_path_choices);

  auto [max_gangsters_index, max_gangsters_length, increasing_order] =
      FindLongestSequenc(greater_numbers, lesser_numbers);

  auto result_gangsters = RestoreBestSequence(
      max_gangsters_index, increasing_order, greater_path_choices,
      lesser_path_choices, numbers_of_gangsters);

  std::cout << max_gangsters_length << std::endl;
  for (const auto& value : result_gangsters) {
    std::cout << value << " ";
  }

  std::cout << std::endl;
}

void CalculateLongestSequences(const std::vector<int>& gangsters,
                               std::vector<int>& greater_numbers,
                               std::vector<int>& greater_path_choices,
                               std::vector<int>& lesser_numbers,
                               std::vector<int>& lesser_path_choices) {
  for (std::size_t right = 0; right < gangsters.size(); ++right) {
    for (std::size_t left = 0; left < right; ++left) {
      if (gangsters[right] > gangsters[left] &&
          greater_numbers[right] < lesser_numbers[left] + 1) {
        greater_numbers[right] = lesser_numbers[left] + 1;
        greater_path_choices[right] = left;
      }
      if (gangsters[right] < gangsters[left] &&
          lesser_numbers[right] < greater_numbers[left] + 1) {
        lesser_numbers[right] = greater_numbers[left] + 1;
        lesser_path_choices[right] = left;
      }
    }
  }
}

std::tuple<int, int, bool> FindLongestSequenc(
    const std::vector<int>& greater_numbers,
    const std::vector<int>& lesser_numbers) {
  int max_gangsters_index = 0;
  int max_gangsters_length = 0;
  bool increasing_order = false;
  for (int i = 0; i < greater_numbers.size(); ++i) {
    if (greater_numbers[i] > max_gangsters_length) {
      max_gangsters_length = greater_numbers[i];
      max_gangsters_index = i;
      increasing_order = true;
    }
    if (lesser_numbers[i] > max_gangsters_length) {
      max_gangsters_length = lesser_numbers[i];
      max_gangsters_index = i;
      increasing_order = false;
    }
  }

  return {max_gangsters_index, max_gangsters_length, increasing_order};
}

std::vector<int> RestoreBestSequence(int start_index, bool increasing_order,
                                     std::vector<int>& greater_path_choices,
                                     std::vector<int>& lesser_path_choices,
                                     const std::vector<int>& gangsters) {
  std::vector<int> result_gangsters;
  int current_index = start_index;
  while (current_index != -1) {
    result_gangsters.push_back(gangsters[current_index]);
    current_index = increasing_order ? greater_path_choices[current_index]
                                     : lesser_path_choices[current_index];
    increasing_order = !increasing_order;
  }
  std::reverse(result_gangsters.begin(), result_gangsters.end());

  return result_gangsters;
}
