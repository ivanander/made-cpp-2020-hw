#pragma once

#include <cmath>
#include <algorithm>
#include <iostream>
#include <optional>
#include <tuple>
#include <vector>

inline double sqr(double val) { return val * val; }

bool are_equals(double a, double b, double eps = 1e-7) {
  return abs(a - b) <= eps;
}

struct Point {
  Point(double x, double y) : x(x), y(y) {}

  double x = 0;
  double y = 0;

  bool operator==(const Point& other) const {
    return are_equals(x, other.x) && are_equals(y, other.y);
  }

  bool operator!=(const Point& other) const {
    return !(*this == other);
  }

  double distance(const Point& other) const {
    return sqrt(sqr(x - other.x) + sqr(y - other.y));
  };
};

struct Line {
  Line(const Point& p1, const Point& p2)
    : coef_((p2.y - p1.y) / (p2.x - p1.x)), shift_(p1.y - coef_ * p1.x) {}

  Line(double coef, double shift)
    : coef_(coef), shift_(shift) {}

  Line(const Point& p1, double coef)
    : coef_(coef), shift_(p1.y - coef_ * p1.x)  {}

  bool operator==(const Line& other) const {
    return are_equals(coef_, other.coef_) && are_equals(shift_, other.shift_);
  }

  bool operator!=(const Line& other) const {
    return !(*this == other);
  }

  const double coef_, shift_;
};

class Shape {
 public:
  Shape(std::vector<Point> points) : points_(std::move(points)) {}
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool operator==(const Shape& other) const = 0;
  virtual bool operator!=(const Shape& other) const = 0;
  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflex(const Point& center) = 0;
  virtual void reflex(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;

 protected:
  std::vector<Point> points_;
};

class Polygon : public Shape {
 public:
  Polygon(std::vector<Point> vertices) : Shape(std::move(vertices)) {}

  double perimeter() const override {
    double res = 0;
    for (int i = 0; i < points_.size(); ++i) {
      res += points_[i].distance(points_[(i + 1) % points_.size()]);
    }
    return res;
  }

  double area() const override {
    double s = 0.0;
    for (int i = 0, j = points_.size() - 1; i < points_.size(); ++i, j = i) {
      s += (points_[j].x + points_[i].x) * (points_[j].y - points_[i].y);
    }
    return abs(s / 2.0);
  }
  
  bool operator==(const Shape& other) const {
    Polygon figure({});

    try {
      figure = dynamic_cast<const Polygon&>(other);
    } catch(const std::exception& e) {
      return false;
    }

    if (points_.size() == figure.points_.size()) {
      if (points_.empty()) {
        return true;
      }
      std::vector<Point> duplicated(points_);
      duplicated.insert(duplicated.end(), duplicated.begin(), duplicated.end());

      auto it1 = std::search(duplicated.begin(), duplicated.end(), figure.points_.begin(), figure.points_.end());
      if (it1 != duplicated.end()) {
        return true;
      }

      auto it2 = std::search(duplicated.rbegin(), duplicated.rend(), figure.points_.begin(), figure.points_.end());
      if (it2 != duplicated.rend()) {
        return true;
      }
    }
    return false;
  }


  bool operator!=(const Shape& other) const override {
    return !(*this == other);
  }

  void rotate(const Point& center, double angle) override {
    angle *= M_PI / 180;
    for (auto& point : points_) {
      Point radius_vec(point.x - center.x, point.y - center.y);
      point = { center.x + radius_vec.x * cos(angle) - radius_vec.y * sin(angle),
                center.y + radius_vec.x * sin(angle) + radius_vec.y * cos(angle) };
    }
  }

  void reflex(const Point& center) override {

  }

  void reflex(const Line& axis) override {
    for (Point& pt : points_) {
      double coef = (pt.x + (pt.y - axis.shift_) * axis.coef_) / (1 + sqr(axis.coef_));
      pt = { 2 * coef - pt.x,
             2 * coef * axis.coef_ - pt.y + 2 * axis.shift_ };
    }
  }

  void scale(const Point& center, double coefficient) override {
    for (auto& point : points_) {
      Point radius_vec(point.x - center.x, point.y - center.y);
      point = { center.x + radius_vec.x * coefficient,
                center.y + radius_vec.y * coefficient };
    }
  }

  size_t verticesCount() const {
    return points_.size();
  }

  std::vector<Point> getVertices() const {
    return points_;
  }

 protected:
  std::vector<Point> points_;
};

class Ellipse : public Polygon {
 public:
  Ellipse(const Point& p1, const Point& p2, double distance)
      : Polygon({p1, p2}),
        center_({(p1.x + p2.x) / 2, (p1.y + p2.y) / 2 }),
        semifoc(distance / 2) {}

  double perimeter() const override {
    double e = eccentricity();
    double a = semifoc / e;
    double b = a * sqrt(1 - sqr(e));
    return 2.0 * M_PI * sqrt((sqr(a) + sqr(b)) / 2.0);
  }

  double area() const override {
    double e = eccentricity();
    double b = semifoc * sqrt(1 - sqr(e));
    return M_PI * semifoc * b;
  }

  std::pair<Point, Point> focuses() const {
    return { points_.at(0), points_.at(1)};
  }

  Point center() const {
    return center_;
  }

  double eccentricity() const {
    const auto& p1 = points_[0];
    const auto& p2 = points_[1];
    double dist = 0.5 * sqrt(sqr(p1.x - p2.x) +
        sqr(p1.y - p2.y));
    return dist / semifoc;
  }

 protected:
  Point center_;
  double semifoc;
};

class Circle : public Ellipse {
 public:
  Circle(const Point& center, double r)
    : Ellipse(center, center, 2 * r) {}

  double perimeter() const override {
    return 2 * M_PI * semifoc;
  }

  double radius() const {
    return semifoc;
  }
};

class Triangle : public Polygon {
 public:
  Triangle(const Point& p1, const Point& p2, const Point& p3)
    : Polygon({p1, p2, p3}) {}

  Point centroid() const {

  }

  Circle circumscribedCircle() const {
    Point center1 = { (points_[0].x + points_[1].x) / 2,
                      (points_[0].y + points_[1].y) / 2 };

    double coef1 = (points_[0].y - points_[1].y) /
        (points_[0].x - points_[1].x);

    coef1 = -1 / coef1;
    double b01 = center1.y - coef1 * center1.x;

    Point center2 = { (points_[0].x + points_[2].x) / 2,
                      (points_[0].y + points_[2].y) / 2 };

    double coef2 = (points_[0].y - points_[2].y) /
        (points_[0].x - points_[2].x);

    coef2 = -1 / coef2;

    double b02 = center2.y - coef2 * center2.x;
    double x = (b02 - b01) / (coef1 - coef2);
    double y = coef1 * x + b01;
    double R = sqrt(sqr(points_[0].x - x) + sqr(points_[0].y - y));

    return Circle({x, y}, R);
  }

  Circle inscribedCircle() const {
    Point center1 = { (points_[0].x + points_[1].x) / 2,
                      (points_[0].y + points_[1].y) / 2 };
    double coef1 = (points_[0].y - points_[1].y) /
        (points_[0].x - points_[1].x);
    coef1 = -1 / coef1;
    double b01 = center1.y - coef1 * center1.x;

    Point center2 = { (points_[0].x + points_[2].x) / 2,
                      (points_[0].y + points_[2].y) / 2 };
    double coef2 = (points_[0].y - points_[2].y) /
        (points_[0].x - points_[2].x);
    coef2 = -1 / coef2;

    double b02 = center2.y - coef2 * center2.x;
    double x = (b02 - b01) / (coef1 - coef2);
    double y = coef1 * x + b01;
    double R = sqrt(sqr(points_[0].x - x) + sqr(points_[0].y - y));

    return Circle({x, y}, R);
  }

  bool isInside(const Point& p) const {
    const Point& p1 = this->points_[0];
    const Point& p2 = this->points_[1];
    const Point& p3 = this->points_[2];
    double A = ((p2.y - p3.y) * (p.x - p3.x) + (p3.x - p2.x) * (p.y - p3.y)) /
        ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y));
    double B = ((p3.y - p1.y) * (p.x - p3.x) + (p1.x - p3.x) * (p.y - p3.y)) /
        ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y));
    double C = 1.0 - A - B;
    return (A > 0) && (B > 0) && (C > 0);
  }

  Point orthocenter() const {}
  Line EulerLine() const {}
  Circle ninePointsCircle() const {}
};

class Rectangle : public Polygon {
 public:
  Rectangle(const Point& p1, const Point& p2, double k) : Polygon({}) {
    double dist = p2.distance(p1);
    Point diff((p2.x - p1.x) / dist, (p2.y - p1.y) / dist);
    k = (k < 1) ? 1 / k : k;
    double phi = atan(k);
    Point vec = { cos(phi) * diff.x - sin(phi) * diff.y,
                  sin(phi) * diff.x + cos(phi) * diff.y };
    points_ = { p1, {p1.x + vec.x, p1.y + vec.y},
                     p2, {p2.x - vec.x, p2.y - vec.y} };
  }

  std::pair<Line, Line> diagonals() const {
    return { {points_[0], points_[2]}, {points_[1], points_[3]} };
  }

  Point center() const {
    return {(points_[0].x + points_[2].x) / 2, (points_[0].y + points_[2].y) / 2 };
  }
 private:
};

class Square : public Rectangle {
 public:
  Square(const Point& p1, const Point& p2) : Rectangle({p1, p2, 1}) {}

  Circle circumscribedCircle() const {
    Point center = Rectangle::center();
    double diag = points_[0].distance(this->points_[2]);
    return Circle(center, diag / 2);
  }

  Circle inscribedCircle() const {
    Point center = Rectangle::center();
    double side = points_[0].distance(points_[1]);
    return Circle(center, side / 2);
  }
};