#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <istream>
#include <ranges>
#include <vector>
#include <numbers>
#include <iomanip>

namespace geometry {

class Vector2D {
  friend std::istream& operator>>(std::istream& istream, Vector2D& vector) {
    istream >> vector.x_coordinate_ >> vector.y_coordinate_;
    return istream;
  }

 public:
  Vector2D() = default;

  Vector2D(std::int64_t x_coordinate, std::int64_t y_coordinate)
      : x_coordinate_(x_coordinate), y_coordinate_(y_coordinate) {}

  Vector2D& operator+=(const Vector2D& other) {
    x_coordinate_ += other.x_coordinate_;
    y_coordinate_ += other.y_coordinate_;
    return *this;
  }

  Vector2D operator+(const Vector2D& other) const {
    Vector2D result(*this);
    result += other;
    return result;
  }

  Vector2D& operator-=(const Vector2D& other) {
    x_coordinate_ -= other.x_coordinate_;
    y_coordinate_ -= other.y_coordinate_;
    return *this;
  }

  Vector2D operator-(const Vector2D& other) const {
    Vector2D result(*this);
    result -= other;
    return result;
  }

  Vector2D& operator*=(std::int64_t scalar) {
    x_coordinate_ *= scalar;
    y_coordinate_ *= scalar;
    return *this;
  }

  Vector2D operator*(std::int64_t scalar) const {
    Vector2D result(*this);
    result *= scalar;
    return result;
  }

  std::int64_t Dot(const Vector2D& other) const {
    return x_coordinate_ * other.x_coordinate_ +
           y_coordinate_ * other.y_coordinate_;
  }

  std::int64_t Cross(const Vector2D& other) const {
    return x_coordinate_ * other.y_coordinate_ -
           y_coordinate_ * other.x_coordinate_;
  }

  std::int64_t LengthSquared() const { return Dot(*this); }

  std::int64_t GetX() const { return x_coordinate_; }
  std::int64_t GetY() const { return y_coordinate_; }

  double PolarAngle() const {
    double angle = std::atan2(y_coordinate_, x_coordinate_);
    if (angle < 0.0) {
      angle += 2.0 * std::numbers::pi;
    }
    return angle;
  }

 private:
  std::int64_t x_coordinate_{0};
  std::int64_t y_coordinate_{0};
};

class Point2D {
  friend std::istream& operator>>(std::istream& istream, Point2D& point) {
    istream >> point.x_coordinate_ >> point.y_coordinate_;
    return istream;
  }

 public:
  Point2D() = default;

  explicit Point2D(const Vector2D& position)
      : x_coordinate_(position.GetX()), y_coordinate_(position.GetY()) {}

  Point2D(std::int64_t x_coordinate, std::int64_t y_coordinate)
      : x_coordinate_(x_coordinate), y_coordinate_(y_coordinate) {}

  Point2D operator+(const Vector2D& vector) const {
    return Point2D(GetPosition() + vector);
  }

  Vector2D operator-(const Point2D& other) const {
    return GetPosition() - other.GetPosition();
  }

  Point2D operator-() const {
    return Point2D(-x_coordinate_, -y_coordinate_);
  }

  Vector2D GetPosition() const {
    return Vector2D(x_coordinate_, y_coordinate_);
  }

  double Distance(const Point2D& other) const {
    return std::hypot(GetX() - other.GetX(), GetY() - other.GetY());
  }

  std::int64_t GetX() const { return x_coordinate_; }
  std::int64_t GetY() const { return y_coordinate_; }

  auto operator<=>(const Point2D& other) const {
    if (GetX() != other.GetX()) {
      return GetX() <=> other.GetX();
    }
    return GetY() <=> other.GetY();
  }

  bool operator==(const Point2D& other) const {
    return GetX() == other.GetX() && GetY() == other.GetY();
  }

 private:
  std::int64_t x_coordinate_{0};
  std::int64_t y_coordinate_{0};
};

class Segment {
 public:
  Segment() = default;

  Segment(Point2D start, Point2D end) : start_(start), end_(end) {}

  Vector2D ToVector() const { return end_ - start_; }

  double Distance(const Point2D& target_point) const {
    Vector2D segment_vector = ToVector();
    Vector2D vector_to_target = target_point - start_;

    std::int64_t segment_length_squared = segment_vector.LengthSquared();

    if (segment_length_squared == 0) {
      return target_point.Distance(start_);
    }

    double projection_coeff =
        vector_to_target.Dot(segment_vector) / static_cast<double>(segment_length_squared);
    projection_coeff = std::max(0.0, std::min(1.0, projection_coeff));

    double closest_x = start_.GetX() + segment_vector.GetX() * projection_coeff;
    double closest_y = start_.GetY() + segment_vector.GetY() * projection_coeff;

    return std::hypot(target_point.GetX() - closest_x, target_point.GetY() - closest_y);
  }

 private:
  Point2D start_;
  Point2D end_;
};

class Polygon {
 public:
  Polygon() = default;

  explicit Polygon(std::size_t size) 
      : polygon_(size) {}

  explicit Polygon(std::vector<Point2D> polygon)
      : polygon_(std::move(polygon)) {}

  auto begin() const { return polygon_.begin(); }
  auto end() const { return polygon_.end(); }

  std::size_t size() const { return polygon_.size(); }

  const Point2D& operator[](std::size_t index) const { return polygon_[index]; }

  Polygon operator-() const {
    std::vector<Point2D> negated;
    negated.reserve(size());
    for (const auto& point : polygon_) {
      negated.push_back(-point);
    }
    return Polygon(std::move(negated));
  }

  void NormalizeBottomLeft() {
    auto min_iterator = std::ranges::min_element(polygon_, [](const Point2D& point_lhs, const Point2D& point_rhs) {
      if (point_lhs.GetY() != point_rhs.GetY()) {
        return point_lhs.GetY() < point_rhs.GetY();
      }
      return point_lhs.GetX() < point_rhs.GetX();
    });
    std::ranges::rotate(polygon_, min_iterator);
  }

  std::int64_t DoubleArea() const {
    size_t points_count = size();
    if (points_count < 3) { 
      return 0;
    }

    std::int64_t area_double = 0;
    for (std::size_t index : std::views::iota(0uz, points_count)) {
      std::size_t next_index = (index + 1) % points_count;
      area_double += ((*this)[index].GetX() * (*this)[next_index].GetY()) - 
                     ((*this)[next_index].GetX() * (*this)[index].GetY());
    }
    return std::abs(area_double);
  }

 private:
  std::vector<Point2D> polygon_;
};

class MinkowskiSum {
 public:
  MinkowskiSum(Polygon lhs, Polygon rhs)
      : lhs_(std::move(lhs)), rhs_(std::move(rhs)) {
    lhs_.NormalizeBottomLeft();
    rhs_.NormalizeBottomLeft();
  }

  Polygon Build() const {
    auto merged_edges = MergeEdges(ExtractEdges(lhs_), ExtractEdges(rhs_));

    std::vector<Point2D> result;
    result.reserve(merged_edges.size());

    Point2D current_point = lhs_[0] + rhs_[0].GetPosition();
    result.push_back(current_point);

    std::ranges::for_each(
        merged_edges | std::views::take(merged_edges.size() - 1),
        [&current_point, &result](const Vector2D& edge_vector) {
            current_point = current_point + edge_vector;
            result.push_back(current_point);
        }
    );

    return Polygon(std::move(result));
  }

 private:
  static std::vector<Vector2D> ExtractEdges(const Polygon& polygon) {
    std::vector<Vector2D> edges;
    edges.reserve(polygon.size());
    for (std::size_t index : std::views::iota(0uz, polygon.size())) {
      edges.push_back(polygon[(index + 1) % polygon.size()] - polygon[index]);
    }
    return edges;
  }

  static std::vector<Vector2D> MergeEdges(const std::vector<Vector2D>& lhs_edges, 
                                          const std::vector<Vector2D>& rhs_edges) {
    std::vector<Vector2D> merged;
    merged.reserve(lhs_edges.size() + rhs_edges.size());
    
    std::ranges::merge(
        lhs_edges, rhs_edges, 
        std::back_inserter(merged),
        [](const Vector2D& edge_vector_lhs, const Vector2D& edge_vector_rhs) {
            return edge_vector_lhs.PolarAngle() < edge_vector_rhs.PolarAngle();
        }
    );
    
    return merged;
  }

  Polygon lhs_;
  Polygon rhs_;
};

std::istream& operator>>(std::istream& istream, Polygon& polygon) {
  std::vector<geometry::Point2D> points(polygon.size());
  for (auto& point : points) {
    std::cin >> point;
  }

  polygon = geometry::Polygon(std::move(points));

  return istream;
}

}  // namespace geometry

double CalculateMinimumTime(const geometry::Polygon& difference) {
  double min_distance = 1e12; 
  geometry::Point2D origin_point(0, 0);

  std::size_t diff_vertex_count = difference.size();
  for (std::size_t index : std::views::iota(0uz, diff_vertex_count)) {
    std::size_t next_index = (index + 1) % diff_vertex_count;
    geometry::Segment border_edge(difference[index], difference[next_index]);
    
    min_distance = std::min(min_distance, border_edge.Distance(origin_point));
  }

  return std::max(0.0, min_distance - 60.0);
}

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::size_t airport_vertex_count;
  std::size_t cloud_vertex_count;

  std::cin >> airport_vertex_count >> cloud_vertex_count;

  geometry::Polygon airport(airport_vertex_count);
  std::cin >> airport;

  geometry::Polygon cloud(cloud_vertex_count);
  std::cin >> cloud;

  geometry::Polygon minkowski_diff = geometry::MinkowskiSum(airport, -cloud).Build();
  auto minimum_time = CalculateMinimumTime(minkowski_diff);
  std::cout << std::fixed << std::setprecision(8) << minimum_time << "\n";

  return 0;
}
