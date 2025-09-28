#include <algorithm>
#include <iostream>
#include <vector>

int main() {
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

  for (int right = 0; right < number_of_districts; ++right) {
    for (int left = 0; left < right; ++left) {
      if (numbers_of_gangsters[right] > numbers_of_gangsters[left] &&
          greater_numbers[right] < lesser_numbers[left] + 1) {
        greater_numbers[right] = lesser_numbers[left] + 1;
        greater_path_choices[right] = left;
      }
      if (numbers_of_gangsters[right] < numbers_of_gangsters[left] &&
          lesser_numbers[right] < greater_numbers[left] + 1) {
        lesser_numbers[right] = greater_numbers[left] + 1;
        lesser_path_choices[right] = left;
      }
    }
  }

  int max_gangsters_index = 0;
  int max_gangsters_length = 0;
  bool increasing_order = false;
  for (int i = 0; i < number_of_districts; ++i) {
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

  std::vector<int> result_gangsters;
  int current_index = max_gangsters_index;
  while (current_index != -1) {
    result_gangsters.push_back(numbers_of_gangsters[current_index]);
    current_index = increasing_order ? greater_path_choices[current_index]
                                     : lesser_path_choices[current_index];
    increasing_order = !increasing_order;
  }
  std::reverse(result_gangsters.begin(), result_gangsters.end());

  std::cout << max_gangsters_length << std::endl;
  for (const auto& value : result_gangsters) {
    std::cout << value << " ";
  }

  std::cout << std::endl;
}

