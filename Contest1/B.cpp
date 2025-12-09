#include <array>
#include <cstdint>
#include <initializer_list>
#include <iostream>

constexpr uint64_t kLengthOfJump = 5;
constexpr uint64_t kModConst = 1000003;

template <typename Type, size_t row_count, size_t column_count>
class Matrix {
 public:
  static constexpr size_t kRows = row_count;
  static constexpr size_t kColumns = column_count;

  constexpr Matrix() : elements_() {}

  constexpr Matrix(
      std::initializer_list<std::initializer_list<Type>> elements_init)
      : elements_() {
    size_t row_index = 0;
    for (const auto& row : elements_init) {
      size_t column_index = 0;
      for (const auto& element : row) {
        elements_[row_index][column_index] = element % kModConst;
        column_index++;
      }
      row_index++;
    }
  }

  template <size_t other_column_count>
  constexpr Matrix<Type, row_count, other_column_count> operator*(
      const Matrix<Type, kColumns, other_column_count>& other) const {
    Matrix<Type, row_count, other_column_count> result;
    for (size_t i = 0; i < kRows; ++i) {
      for (size_t j = 0; j < other.kColumns; ++j) {
        Type result_sum = 0;
        for (size_t k = 0; k < kColumns; ++k) {
          result_sum =
              (result_sum + GetElement(i, k) * other.GetElement(k, j)) %
              kModConst;
        }
        result.GetElement(i, j) = result_sum;
      }
    }
    return result;
  }

  constexpr Matrix<Type, kColumns, kRows> Transpose() const {
    Matrix<Type, kColumns, kRows> result;
    for (size_t i = 0; i < kRows; ++i) {
      for (size_t j = 0; j < kColumns; ++j) {
        result.GetElement(j, i) = GetElement(i, j);
      }
    }
    return result;
  }

  constexpr Matrix Identity() const {
    Matrix<Type, kRows, kColumns> result;
    for (size_t i = 0; i < kRows; ++i) {
      result.GetElement(i, i) = 1;
    }
    return result;
  }

  constexpr Matrix Pow(uint64_t power) const {
    if (power == 0) {
      return Identity();
    }

    if (power == 1) {
      return *this;
    }

    if (power % 2 == 0) {
      Matrix half_power_matrix = Pow(power / 2);
      return half_power_matrix * half_power_matrix;
    }

    Matrix minus_one_power_matrux = Pow(power - 1);
    return minus_one_power_matrux * (*this);
  }

  constexpr Type& GetElement(size_t row, size_t column) {
    return elements_[row][column];
  }

  constexpr const Type& GetElement(size_t row, size_t column) const {
    return elements_[row][column];
  }

  constexpr void Print() const {
    for (size_t i = 0; i < row_count; ++i) {
      for (size_t j = 0; j < column_count; ++j) {
        std::cout << GetElement(i, j) << " ";
      }
      std::cout << std::endl;
    }
  }

 private:
  std::array<std::array<Type, column_count>, row_count> elements_;
};

constexpr Matrix<uint64_t, 1, 5> kFiveAnswers = {{8, 4, 2, 1, 1}};
constexpr Matrix<uint64_t, 5, 5> kGrassHopperMatrix = {{1, 1, 1, 1, 1},
                                                       {1, 0, 0, 0, 0},
                                                       {0, 1, 0, 0, 0},
                                                       {0, 0, 1, 0, 0},
                                                       {0, 0, 0, 1, 0}};

int main() {
  uint64_t number = 0;
  std::cin >> number;

  if (number <= kLengthOfJump) {
    std::cout << kFiveAnswers.GetElement(0, kLengthOfJump - number)
              << std::endl;
    return 0;
  }

  Matrix result = kGrassHopperMatrix.Pow(number - 5) * kFiveAnswers.Transpose();
  std::cout << result.GetElement(0, 0) << std::endl;
}

