#include <algorithm>
#include <iostream>
#include <vector>

class ResourceBackpack {
 public:
  ResourceBackpack(int max_resources, std::vector<int> resources,
                   std::vector<int> rewards)
      : max_resources_(max_resources),
        mission_resources_(resources),
        mission_rewards_(rewards) {
    std::vector<std::vector<int>> temp_backpack(
        resources.size() + 1, std::vector<int>(max_resources + 1, 0));
    resources_backpack_.swap(temp_backpack);
  }

  void CalculateBackpack() {
    for (std::size_t i = 1; i <= mission_resources_.size(); ++i) {
      for (int j = 0; j <= max_resources_; ++j) {
        resources_backpack_[i][j] = resources_backpack_[i - 1][j];
        if (mission_resources_[i - 1] <= j) {
          resources_backpack_[i][j] = std::max(
              resources_backpack_[i][j],
              resources_backpack_[i - 1][j - mission_resources_[i - 1]] +
                  mission_rewards_[i - 1]);
        }
      }
    }
  }

  std::vector<std::size_t> RestorePath() {
    std::vector<std::size_t> result_path;
    std::size_t current_resource = mission_resources_.size();
    int current_mass = max_resources_;
    while (current_resource > 0 && current_mass > 0) {
      if (resources_backpack_[current_resource][current_mass] !=
          resources_backpack_[current_resource - 1][current_mass]) {
        result_path.push_back(current_resource);
        current_mass -= mission_resources_[current_resource - 1];
      }
      current_resource--;
    }
    std::reverse(result_path.begin(), result_path.end());

    return result_path;
  }

  void PrintBackpack() {
    for (const auto& row : resources_backpack_) {
      for (const auto& value : row) {
        std::cout << value << "\t";
      }
      std::cout << "\n";
    }
  }

 private:
  int max_resources_;

  std::vector<int> mission_resources_;
  std::vector<int> mission_rewards_;

  std::vector<std::vector<int>> resources_backpack_;
};

int main() {
  int missions_number;
  int max_resources;
  std::cin >> missions_number >> max_resources;

  std::vector<int> mission_resources(missions_number);
  for (int i = 0; i < missions_number; ++i) {
    std::cin >> mission_resources[i];
  }

  std::vector<int> mission_rewards(missions_number);
  for (int i = 0; i < missions_number; ++i) {
    std::cin >> mission_rewards[i];
  }

  ResourceBackpack resources_backpack(max_resources, mission_resources,
                                      mission_rewards);
  resources_backpack.CalculateBackpack();
  for (const auto& value : resources_backpack.RestorePath()) {
    std::cout << value << "\n";
  }

  return 0;
}

