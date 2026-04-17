#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <bit>

constexpr std::int32_t kMod = 7'340'033;
constexpr std::int32_t kFirstOrderRoot = 3;

struct Zp {
  std::int32_t value;

  Zp(std::int64_t val = 0) 
  : value(val % kMod) {
    if (value < 0) {
      value += kMod;
    }
  }

  Zp operator+(const Zp& other) {
    return Zp(value + other.value);
  }

  Zp operator-(const Zp& other) {
    return Zp(value - other.value);
  }

  Zp operator*(const Zp& other) {
    return Zp(1LL * value * other.value);
  }

  Zp Pow(std::int64_t power) const {
    Zp result = 1;
    Zp self = *this; 

    for (; power > 0; power >>= 1, self = self * self) {
      if (power % 2 == 1) {
        result = result * self;
      }
    }

    return result;
  }

  Zp Inverse() const {
    return Pow(kMod - 2);
  }
};

struct VectorZp : std::vector<Zp> {
  friend VectorZp operator-(VectorZp lhs, const VectorZp& rhs) {
    return lhs -= rhs;
  }

  using std::vector<Zp>::vector;

  VectorZp(Zp value) : std::vector<Zp>({value}) {}

  VectorZp Slice(std::size_t index) const {
    return VectorZp(begin(), begin() + std::min(size(), index));
  }

  VectorZp operator-() const {
    VectorZp result = *this;
    for (auto& element : result) {
      element = Zp(0) - element;
    }
    return result;
  }

  VectorZp& operator-=(const VectorZp& other) {
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
static std::vector<Zp> root_cache(kMaxLog + 1);
static std::vector<Zp> inverse_root_cache(kMaxLog + 1);

static const Zp kInverseTwo = Zp(2).Inverse();

struct NumberTheoreticTransform {
  static void Prepare() {
    for (int i = 0; i <= kMaxLog; ++i) {
      root_cache[i] = Zp(kFirstOrderRoot).Pow((kMod - 1) >> i);
      inverse_root_cache[i] = root_cache[i].Inverse();
    }
  }

  static void Transform(VectorZp& coefficients, bool invert) {
    std::size_t coeff_size = coefficients.size();
    if (coeff_size <= 1) {
      return;
    }

    VectorZp even_coefficients(coeff_size / 2);
    VectorZp odd_coefficients(coeff_size / 2);

    for (std::size_t i = 0; i < coeff_size / 2; ++i) {
      even_coefficients[i] = coefficients[i * 2];
      odd_coefficients[i] = coefficients[(i * 2) + 1];
    }

    Transform(even_coefficients, invert);
    Transform(odd_coefficients, invert);

    std::size_t log_coeff_size = std::__countr_zero(coeff_size);

    Zp current_root = 1;
    Zp step_root = invert ? inverse_root_cache[log_coeff_size] : root_cache[log_coeff_size];

    for (std::size_t i = 0; i < coeff_size / 2; ++i) {
      Zp turned_odd_value = current_root * odd_coefficients[i];

      coefficients[i] = even_coefficients[i] + turned_odd_value;
      coefficients[i + (coeff_size / 2)] = even_coefficients[i] - turned_odd_value;

      if (invert) {
        coefficients[i] = coefficients[i] * kInverseTwo;
        coefficients[i + (coeff_size / 2)] = coefficients[i + (coeff_size / 2)] * kInverseTwo;
      }

      current_root = current_root * step_root;
    }

  }

};

class Polynomial {
public:
  explicit Polynomial(const VectorZp& coefficients)
  : coefficients_(coefficients) {}

  static Polynomial Multiply(Polynomial lhs, Polynomial rhs, std::size_t limit) {
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

  Polynomial Inverse(std::size_t limit) const {
    if (coefficients_.empty() 
      || coefficients_[0].value == 0) {
      throw std::runtime_error("The ears of a dead donkey");
    }

    Polynomial Q(VectorZp(coefficients_[0].Inverse()));
    for (std::size_t i = 2; i <= (limit << 1); i <<= 1) {
      Polynomial P_i(coefficients_.Slice(i));
      Q = Multiply(Q, Polynomial(Zp(2) - Multiply(P_i, Q, i).coefficients_), i);
    }

    Q.Resize(limit);
    return Q;
  }

  Zp& operator[](std::size_t index) {
    return coefficients_[index];
  } 

  std::size_t Degree() const {
    return coefficients_.size();
  }

  void Resize(std::size_t size) {
    coefficients_.resize(size);
  }

  void Slice(std::size_t limit) {
    coefficients_.resize(limit, Zp(0));
  }

private:
  VectorZp coefficients_;
};

int main() {
  NumberTheoreticTransform::Prepare();

  std::size_t lhs_degree; 
  std::size_t rhs_degree;
  std::cin >> lhs_degree >> rhs_degree;

  VectorZp polynom_coefficients(lhs_degree + 1, Zp(0));
  for (std::size_t i = 0; i <= rhs_degree; ++i) {
    std::int64_t value = 0;
    std::cin >> value;
    if (i < lhs_degree) {
      polynom_coefficients[i] = Zp(value);
    }
  }

  Polynomial P(polynom_coefficients);
  try {
    Polynomial Q = P.Inverse(lhs_degree);
    for (std::size_t i = 0; i < lhs_degree; ++i) {
      std::cout << Q[i].value << " ";
    }
    std::cout << "\n";
  } catch (std::runtime_error error) {
    std::cout << error.what() <<  "\n";
  }
}


