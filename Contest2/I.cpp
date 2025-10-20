#include <iostream>
#include <istream>
#include <string>
#include <vector>

class SmallBigInt {
  friend std::istream& operator>>(std::istream& istream, SmallBigInt& number) {
    istream >> number.value_;
    if (number.value_.empty()) {
      number.value_ = "0";
    }
    return istream;
  }

 private:
  static const int kBase = 10;

 public:
  SmallBigInt() = default;

  SmallBigInt(const std::string& str) : value_(str.empty() ? "0" : str) {}

  bool IsZero() const { return value_ == "0"; }

  bool IsOdd() const { return GetDigit(value_.back()) % 2 != 0; }

  SmallBigInt SubstractOne() const {
    if (value_ == "0") {
      return SmallBigInt();
    }

    std::string new_value = value_;
    int remainder = -1;
    for (std::size_t i = value_.length() - 1; i >= 0; --i) {
      int digit = GetDigit(new_value[i]) + remainder;
      if (digit < 0) {
        new_value[i] = GetChar(digit + kBase);
        remainder = -1;
      } else {
        new_value[i] = GetChar(digit);
        break;
      }
    }

    if (new_value.length() > 1 && new_value[0] == '0') {
      new_value.erase(0, 1);
    }

    return SmallBigInt(new_value);
  }

  void DivideBy2() {
    if (value_ == "0") {
      return;
    }

    std::string new_value;
    int remainder = 0;

    for (const auto& symbol : value_) {
      int digit = GetDigit(symbol);
      int current_number = (remainder * kBase) + digit;
      int result_digit = current_number / 2;
      remainder = current_number % 2;

      if (!new_value.empty() || result_digit != 0) {
        new_value += std::to_string(result_digit);
      }
    }

    value_ = new_value;
    value_ = new_value.empty() ? "0" : new_value;
  }

 private:
  static int GetDigit(char symbol) { return symbol - '0'; }

  static char GetChar(int digit) { return static_cast<char>(digit) + '0'; }

  std::string value_ = "0";
};

class MatrixWithBigInt {
  friend MatrixWithBigInt FastPow(MatrixWithBigInt base, SmallBigInt degree) {
    MatrixWithBigInt result(base.size_, base.module_);
    result.SetIdentity();
    while (!degree.IsZero()) {
      if (degree.IsOdd()) {
        result = result * base;
      }
      base = base * base;
      degree.DivideBy2();
    }
    return result;
  }

 public:
  MatrixWithBigInt(std::size_t size, long long module)
      : size_(size), module_(module) {
    data_.assign(size, std::vector<long long>(size, 0));
  }

  long long& operator()(std::size_t row, std::size_t column) {
    return data_[row][column];
  }

  const long long& operator()(std::size_t row, std::size_t column) const {
    return data_[row][column];
  }

  void SetIdentity() {
    for (std::size_t i = 0; i < size_; ++i) {
      for (std::size_t j = 0; j < size_; ++j) {
        data_[i][j] = (i == j) ? 1 : 0;
      }
    }
  }

  MatrixWithBigInt operator*(const MatrixWithBigInt& other) const {
    MatrixWithBigInt result(size_, module_);
    for (std::size_t i = 0; i < size_; ++i) {
      for (std::size_t j = 0; j < size_; ++j) {
        long long element = 0;
        for (std::size_t k = 0; k < size_; ++k) {
          element = (element + data_[i][k] * other.data_[k][j]) % module_;
        }
        result.data_[i][j] = element;
      }
    }
    return result;
  }

  long long SumElements() const {
    long long result = 0;
    for (std::size_t i = 0; i < size_; ++i) {
      for (std::size_t j = 0; j < size_; ++j) {
        result = (result + data_[i][j]) % module_;
      }
    }
    return result;
  }

 private:
  std::vector<std::vector<long long>> data_;
  std::size_t size_;
  long long module_;
};

bool CheckTransition(std::size_t previous_mask, std::size_t current_mask,
                     std::size_t row_count);
long long CountVariants(const SmallBigInt& columns, const long long& rows,
                        std::size_t module);

int main() {
  SmallBigInt columns;
  long long rows = 0;
  long long module = 0;

  std::cin >> columns >> rows >> module;
  std::cout << CountVariants(columns, rows, module) << "\n";

  return 0;
}

bool CheckTransition(std::size_t previous_mask, std::size_t current_mask,
                     std::size_t row_count) {
  for (std::size_t row = 0; row < row_count - 1; ++row) {
    std::size_t prev_up = (previous_mask >> row) & 1;
    std::size_t prev_bottom = (previous_mask >> (row + 1)) & 1;
    std::size_t curr_up = (current_mask >> row) & 1;
    std::size_t curr_bottom = (current_mask >> (row + 1)) & 1;

    if (prev_up == prev_bottom && prev_up == curr_up &&
        prev_up == curr_bottom) {
      return false;
    }
  }
  return true;
}

long long CountVariants(const SmallBigInt& columns, const long long& rows,
                        std::size_t module) {
  std::size_t matrix_size = 1 << rows;

  if (columns.SubstractOne().IsZero()) {
    return matrix_size % module;
  }

  SmallBigInt power = columns.SubstractOne();
  MatrixWithBigInt transition_matrix(matrix_size, module);
  for (std::size_t previous_mask = 0; previous_mask < matrix_size;
       ++previous_mask) {
    for (std::size_t current_mask = 0; current_mask < matrix_size;
         ++current_mask) {
      if (CheckTransition(previous_mask, current_mask, rows)) {
        transition_matrix(current_mask, previous_mask) = 1;
      }
    }
  }

  MatrixWithBigInt result_matrix = FastPow(transition_matrix, power);
  return result_matrix.SumElements();
}

