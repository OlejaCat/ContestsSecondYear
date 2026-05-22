#include <bit>
#include <cstdint>
#include <iostream>
#include <istream>
#include <stdexcept>
#include <vector>

constexpr std::int32_t kMod = 7'340'033;
constexpr std::int32_t kFirstOrderRoot = 3;

struct NumberModulo {
  std::int32_t value;

  NumberModulo(std::int64_t val = 0) : value(val % kMod) {
    if (value < 0) {
      value += kMod;
    }
  }

  NumberModulo operator+(const NumberModulo& other) {
    return NumberModulo(value + other.value);
  }

  NumberModulo operator-(const NumberModulo& other) {
    return NumberModulo(value - other.value);
  }

  NumberModulo operator*(const NumberModulo& other) {
    return NumberModulo(1LL * value * other.value);
  }
};

inline NumberModulo Pow(NumberModulo number, std::int64_t power) {
  NumberModulo result = 1;

  for (; power > 0; power >>= 1, number = number * number) {
    if (power % 2 == 1) {
      result = result * number;
    }
  }

  return result;
}

inline NumberModulo Inverse(NumberModulo number) { return Pow(number, kMod - 2); }

struct VectorNumberModulo : std::vector<NumberModulo> {
  friend VectorNumberModulo operator-(VectorNumberModulo lhs,
                                      const VectorNumberModulo& rhs) {
    return lhs -= rhs;
  }

  using std::vector<NumberModulo>::vector;

  VectorNumberModulo(NumberModulo value) : std::vector<NumberModulo>({value}) {}

  VectorNumberModulo Slice(std::size_t index) const {
    return VectorNumberModulo(begin(), begin() + std::min(size(), index));
  }

  VectorNumberModulo operator-() const {
    VectorNumberModulo result = *this;
    for (auto& element : result) {
      element = NumberModulo(0) - element;
    }
    return result;
  }

  VectorNumberModulo& operator-=(const VectorNumberModulo& other) {
    if (size() < other.size()) {
      resize(other.size(), 0);
    }

    for (std::size_t i = 0; i < other.size(); ++i) {
      (*this)[i] = (*this)[i] - other[i];
    }

    return *this;
  }
};

static constexpr int kMaxLog = 20;
static std::vector<NumberModulo> root_cache(kMaxLog + 1);
static std::vector<NumberModulo> inverse_root_cache(kMaxLog + 1);

static const NumberModulo kInverseTwo = Inverse(2);

struct NumberTheoreticTransform {
  static void Prepare() {
    for (int i = 0; i <= kMaxLog; ++i) {
      root_cache[i] = Pow(kFirstOrderRoot, (kMod - 1) >> i);
      inverse_root_cache[i] = Inverse(root_cache[i]);
    }
  }

  static void Transform(VectorNumberModulo& coefficients, bool invert) {
    std::size_t coeff_size = coefficients.size();
    if (coeff_size <= 1) {
      return;
    }

    VectorNumberModulo even_coefficients(coeff_size / 2);
    VectorNumberModulo odd_coefficients(coeff_size / 2);

    for (std::size_t i = 0; i < coeff_size / 2; ++i) {
      even_coefficients[i] = coefficients[i * 2];
      odd_coefficients[i] = coefficients[(i * 2) + 1];
    }

    Transform(even_coefficients, invert);
    Transform(odd_coefficients, invert);

    std::size_t log_coeff_size = std::__countr_zero(coeff_size);

    NumberModulo current_root = 1;
    NumberModulo step_root = invert ? inverse_root_cache[log_coeff_size] :
                                      root_cache[log_coeff_size];

    for (std::size_t i = 0; i < coeff_size / 2; ++i) {
      NumberModulo turned_odd_value = current_root * odd_coefficients[i];

      coefficients[i] = even_coefficients[i] + turned_odd_value;
      coefficients[i + (coeff_size / 2)] =
          even_coefficients[i] - turned_odd_value;

      if (invert) {
        coefficients[i] = coefficients[i] * kInverseTwo;
        coefficients[i + (coeff_size / 2)] =
            coefficients[i + (coeff_size / 2)] * kInverseTwo;
      }

      current_root = current_root * step_root;
    }
  }
};

class Polynomial {
 public:
  explicit Polynomial(const VectorNumberModulo& coefficients)
      : coefficients_(coefficients) {}

  static Polynomial Multiply(Polynomial lhs, Polynomial rhs,
                             std::size_t limit) {
    std::size_t degree = 1;
    while (degree < lhs.Degree() + rhs.Degree()) {
      degree <<= 1;
    }

    lhs.Resize(degree);
    rhs.Resize(degree);

    NumberTheoreticTransform::Transform(lhs.coefficients_, false);
    NumberTheoreticTransform::Transform(rhs.coefficients_, false);

    for (std::size_t i = 0; i < degree; ++i) {
      lhs.coefficients_[i] = lhs.coefficients_[i] * rhs.coefficients_[i];
    }

    NumberTheoreticTransform::Transform(lhs.coefficients_, true);

    lhs.Resize(limit);
    return lhs;
  }

  Polynomial InverseNTT(std::size_t limit) const {
    if (coefficients_.empty() || coefficients_[0].value == 0) {
      throw std::runtime_error("The ears of a dead donkey");
    }

    Polynomial Q(VectorNumberModulo(Inverse(coefficients_[0])));
    for (std::size_t i = 2; i <= (limit << 1); i <<= 1) {
      Polynomial P_i(coefficients_.Slice(i));
      Q = Multiply(
          Q, Polynomial(NumberModulo(2) - Multiply(P_i, Q, i).coefficients_),
          i);
    }

    Q.Resize(limit);
    return Q;
  }

  NumberModulo& operator[](std::size_t index) { return coefficients_[index]; }

  std::size_t Degree() const { return coefficients_.size(); }

  void Resize(std::size_t size) { coefficients_.resize(size); }

  void Slice(std::size_t limit) {
    coefficients_.resize(limit, NumberModulo(0));
  }

 private:
  VectorNumberModulo coefficients_;
};

std::istream& operator>>(std::istream& istream, VectorNumberModulo& modulo_vector) {
  for (std::size_t i = 0; i < modulo_vector.size(); ++i) {
    std::int64_t value = 0;
    std::cin >> value;
    modulo_vector[i] = NumberModulo(value);
  }

  return istream;
}

void PrintInversedPolynom(VectorNumberModulo& polynom_coefficients, std::size_t limit) {
  Polynomial P(polynom_coefficients);
  Polynomial Q = P.InverseNTT(limit);
  for (std::size_t i = 0; i < limit; ++i) {
    std::cout << Q[i].value << " ";
  }
  std::cout << "\n";
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  NumberTheoreticTransform::Prepare();

  std::size_t lhs_degree, rhs_degree;
  std::cin >> lhs_degree >> rhs_degree;

  VectorNumberModulo polynom_coefficients(rhs_degree + 1, NumberModulo(0));
  std::cin >> polynom_coefficients;
  polynom_coefficients.Slice(lhs_degree);

  try {
    PrintInversedPolynom(polynom_coefficients, lhs_degree);
  } catch (std::runtime_error error) {
    std::cout << error.what() << "\n";
  }
}
