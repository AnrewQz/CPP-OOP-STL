#include <cmath>
#include <iostream>
#include <vector>

class Line;

struct Point {
  double x = 0;
  double y = 0;
  Point(double x, double y) : x(x), y(y) {}
  Point() {}
  void rotate(const Point& center, double angle);
  void reflect(const Point& center);
  void reflect(const Line& axis);
  void scale(const Point& center, double coefficient);
  Point projection(const Line& line) const;
};

class Line {
private:
  double a = 0;
  double b = 0;
  double c = 0;

public:
  Line(double a, double b, double c) : a(a), b(b), c(c) {}
  Line(const Point& first, const Point& second);
  Line(double angle, double delta) : a(-angle), b(1), c(-delta) {}
  Line(const Point& point, double angle)
      : a(-angle), b(1), c(-a * point.x - point.y) {}
  Line() {}
  double getA() const { return a; }
  double getB() const { return b; }
  double getC() const { return c; }
};

namespace Consts {
  const double kNull = 0.0000001;
  const double kStraight_angle = 180;

  struct Vector {
    double x = 0;
    double y = 0;
    Vector(const Point& begin, const Point& end)
        : x(end.x - begin.x), y(end.y - begin.y) {}
    Vector(const Point& end) : x(end.x), y(end.y) {}
    Vector(double x, double y) : x(x), y(y) {}
    Vector() = default;
    Vector& operator=(const Vector& other) = default;
    ~Vector() = default;
    double Length() { return sqrt(pow(x, 2) + pow(y, 2)); }
  };

  bool Equals(double fst, double snd) { return std::abs(fst - snd) < Consts::kNull; }

  double Scalar(const Vector& first, const Vector& second) {
    return first.x *second.x + first.y *second.y;
  }

  double Pseudo(const Vector& first, const Vector& second) {
    return first.x *second.y - first.y *second.x;
  }

  double Distance(const Point& first, const Point& second = Point(0, 0)) {
    return sqrt(pow(first.x - second.x, 2) + pow(first.y - second.y, 2));
  }

  double Distance(const Point& point, const Line& line) {
    double numerator =
        std::abs(line.getA() * point.x + line.getB() * point.y + line.getC());
    double denominator = sqrt(pow(line.getA(), 2) + pow(line.getB(), 2));
    return numerator / denominator;
  }

  Line Perpendicular(const Point& point, const Line& line) {
    Line result(-line.getB(), line.getA(),
                line.getB() *point.x - line.getA() *point.y);
    return result;
  }

  Point Middle(const Point& first, const Point& second) {
    double x = first.x + second.x;
    double y = first.y + second.y;
    return Point(x / 2, y / 2);
  }
}

class Shape {
public:
  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflect(const Point& center) = 0;
  virtual void reflect(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool operator==(const Shape& another) const = 0;
  virtual bool isCongruentTo(const Shape& another) const = 0;
  virtual bool isSimilarTo(const Shape& another) const = 0;
  virtual bool containsPoint(const Point& another) const = 0;
  virtual ~Shape() = default;
};

class Polygon : public Shape {

protected:
  std::vector<Point> vertices;
  bool isSimilarToOrdrered(const Shape& another) const;

public:
  Polygon() {}
  Polygon(const std::vector<Point>& args) : vertices(args) {}

  template <typename T, typename... Args> 
  void add(T point, Args... args) {
    vertices.push_back(point);
    add(args...);
  }

  void add() {}

  template <typename... Args>
  Polygon(Args... args) { add(args...); }

  size_t verticesCount();
  const std::vector<Point> getVertices() const;
  bool isConvex();
  void rotate(const Point& center, double angle);
  void reflect(const Point& center);
  void reflect(const Line& axis);
  void scale(const Point& center, double coefficient);
  double perimeter() const;
  double area() const;
  bool operator==(const Shape& another) const;
  bool operator!=(const Shape& another) const;
  bool isCongruentTo(const Shape& another) const;
  bool isSimilarTo(const Shape& another) const;
  bool containsPoint(const Point& another) const;
  ~Polygon() = default;
};

class Ellipse : public Shape {
protected:
  Point first;
  Point second;
  double rad;

  double a() const { return rad / 2; }
  double b() const { return sqrt(pow(a(), 2) - pow(c(), 2)); }
  double c() const { return Consts::Distance(first, second) / 2; }

public:
  Ellipse() {}
  Ellipse(Point fst, Point snd, double rad)
      : first(fst), second(snd), rad(rad) {}
  std::pair<Point, Point> focuses() const { return {first, second}; }
  std::pair<Line, Line> directrices() const;
  double eccentricity() const { return c() / a(); };
  void rotate(const Point& center, double angle);
  void reflect(const Point& center);
  void reflect(const Line& axis);
  void scale(const Point& center, double coefficient);
  double perimeter() const;
  double area() const;
  bool operator==(const Shape& another) const;
  bool operator!=(const Shape& another) const;
  bool isCongruentTo(const Shape& another) const;
  bool isSimilarTo(const Shape& another) const;
  bool containsPoint(const Point& another) const;
  ~Ellipse() = default;
};

class Circle : public Ellipse {
public:
  Circle(const Point& center, double rad) : Ellipse(center, center, 2 * rad) {}
  double radius() const { return rad / 2; }
  Point center() const { return first; }
  virtual ~Circle() = default;
};

class Rectangle : public Polygon {
public:
  Rectangle(const Point& first, const Point& second, double coefficient);
  Point center() const { return Consts::Middle(vertices[0], vertices[2]); }
  std::pair<Line, Line> diagonals() const;
  virtual ~Rectangle() = default;
};

class Square : public Rectangle {
public:
  Square(const Point& first, const Point& third)
      : Rectangle(first, third, 1.0) {}
  double side() const { return Consts::Distance(vertices[0], vertices[1]); }
  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
  ~Square() = default;
};

class Triangle : public Polygon {
public:
  using Polygon::Polygon;
  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
  Point circumcenter() const;
  Point centroid() const;
  Point orthocenter() const;
  Line EulerLine() const;
  Circle ninePointsCircle() const;
  Line bisector(int index) const;
  ~Triangle() = default;
};

bool operator==(const Point& first, const Point& second) {
  return Consts::Equals(first.x, second.x) && Consts::Equals(first.y, second.y);
}

bool operator!=(const Point& first, const Point& second) {
  return !(first == second);
}

void Point::rotate(const Point& center, double angle) {
  angle = (angle / Consts::kStraight_angle) * M_PI;
  double prevX = x;
  x = (x - center.x) *cos(angle) - (y - center.y) *sin(angle);
  y = (prevX - center.x) *sin(angle) + (y - center.y) *cos(angle);
  x += center.x;
  y += center.y;
}

void Point::reflect(const Point& center) {
  x = 2 *center.x - x;
  y = 2 *center.y - y;
}

void Point::reflect(const Line& axis) {
  Point temp = projection(axis);
  reflect(temp);
}

void Point::scale(const Point& center, double coefficient) {
  x = center.x + (x - center.x) *coefficient;
  y = center.y + (y - center.y) *coefficient;
}

Point intersection(const Line& first, const Line& second) {
  double mainDet = first.getA() *second.getB() - first.getB() *second.getA();
  Point result;
  double firstDet = first.getB() *second.getC() - first.getC() *second.getB();
  double secondDet =
      first.getC() *second.getA() - first.getA() *second.getC();
  result.x = firstDet / mainDet;
  result.y = secondDet / mainDet;
  return result;
}

Point Point::projection(const Line& line) const {
  return intersection(line, Consts::Perpendicular(*this, line));
}

Line::Line(const Point& first, const Point& second) {
  if (!Consts::Equals(first.x, second.x)) {
    double k = (first.y - second.y) / (second.x - first.x);
    a = k;
    b = 1;
    c = -(first.y + first.x *k);
  } else {
    double k = (second.x - first.x) / (first.y - second.y);
    a = 1;
    b = k;
    c = -(first.x + k *first.y);
  }
}

bool operator==(const Line& first, const Line& second) {
  if (Consts::Equals(first.getA() *second.getB(),
                     first.getB() *second.getA())) {
    if (Consts::Equals(first.getA() *second.getC(),
                       first.getC() *second.getA())) {
      return true;
    }
  }
  return false;
}

bool operator!=(const Line& first, const Line& second) {
  return !(first == second);
}

size_t Polygon::verticesCount() { return vertices.size(); }
const std::vector<Point> Polygon::getVertices() const { return vertices; }

bool Polygon::isConvex() {
  int first = 0;
  int second = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    Point prev = vertices[(vertices.size() + i - 1) % vertices.size()];
    Point next = vertices[(i + 1) % vertices.size()];
    Consts::Vector fst(prev, vertices[i]);
    Consts::Vector snd(vertices[i], next);
    Consts::Pseudo(fst, snd) > 0 ? ++first : ++second;
  }
  return (first == 0 || second == 0);
}

void Polygon::rotate(const Point& center, double angle) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].rotate(center, angle);
  }
}

void Polygon::reflect(const Point& center) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].reflect(center);
  }
}

void Polygon::reflect(const Line& axis) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].reflect(axis);
  }
}

void Polygon::scale(const Point& center, double coefficient) {
  for (size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].scale(center, coefficient);
  }
}

double Polygon::perimeter() const {
  double result = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    Point now = vertices[i];
    Point next = vertices[(i + 1) % vertices.size()];
    result += Consts::Distance(now, next);
  }
  return result;
}

double Polygon::area() const {
  double result = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    Consts::Vector now(vertices[i]);
    Consts::Vector next(vertices[(i + 1) % vertices.size()]);
    result += Consts::Pseudo(now, next);
  }
  return std::abs(result / 2);
}

bool Polygon::operator==(const Shape& another) const {
  Polygon temp = dynamic_cast<const Polygon& >(another);
  if (vertices.size() != temp.vertices.size()) {
    return false;
  }
  for (size_t step = 0; step < vertices.size(); ++step) {
    bool result = true;
    for (size_t i = 0; i < vertices.size(); ++i) {
      if (vertices[i] != temp.vertices[(i + step) % vertices.size()]) {
        result = false;
        break;
      }
    }
    if (result == true) {
      return true;
    }
    result = true;
    for (size_t i = 0; i < vertices.size(); ++i) {
      if (vertices[i] !=
          temp.vertices[(vertices.size() + step - i) % vertices.size()]) {
        result = false;
        break;
      }
    }
    if (result == true) {
      return true;
    }
  }
  return false;
}

bool Polygon::operator!=(const Shape& another) const {
  return !(*this == another);
}

bool Polygon::isSimilarTo(const Shape& another) const {
  std::vector<Point> reversed;
  for (long long i = vertices.size() - 1; i >= 0; --i) {
    reversed.push_back(vertices[i]);
  }
  Polygon rev(reversed);
  return isSimilarToOrdrered(another) || rev.isSimilarToOrdrered(another);
}

bool Polygon::isSimilarToOrdrered(const Shape& another) const {
  const Polygon* temp = dynamic_cast<const Polygon* >(&another);
  if (temp == nullptr) {
    return false;
  }
  if (vertices.size() != temp->vertices.size()) {
    return false;
  }
  for (size_t step = 0; step < vertices.size(); ++step) {
    bool ans = true;
    size_t sz = vertices.size();
    for (size_t i = 0; i < sz; ++i) {
      int prev = (sz + i - 1) % sz;
      int now = i;
      int next = (i + 1) % sz;
      Consts::Vector firstV(vertices[prev], vertices[now]);
      Consts::Vector secondV(vertices[now], vertices[next]);
      double firstCos = Consts::Scalar(firstV, secondV) /
                        (firstV.Length() *secondV.Length());
      Consts::Vector thirdV(temp->vertices[(now + step) % sz],
                            temp->vertices[(prev + step) % sz]);
      Consts::Vector fourthV(temp->vertices[(next + step) % sz],
                             temp->vertices[(now + step) % sz]);
      double secondCos = Consts::Scalar(thirdV, fourthV) /
                         (thirdV.Length() *fourthV.Length());
      if (!Consts::Equals(firstCos, secondCos)) {
        ans = false;
        break;
      }
    }
    if (ans == true) {
      return true;
    }
  }
  return false;
}

bool Polygon::isCongruentTo(const Shape& another) const {
  return (Consts::Equals(area(), another.area())) && isSimilarTo(another);
}

bool Polygon::containsPoint(const Point& another) const {
  int count = 0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    Point first = vertices[i];
    Point second = vertices[(i + 1) % vertices.size()];
    if (Consts::Equals(Consts::Distance(another, first) +
                           Consts::Distance(another, second),
                       Consts::Distance(first, second))) {
      return true;
    }
    if (((first.y - another.y) > 0 && (second.y - another.y) <= 0) ||
        ((first.y - another.y) <= 0 && (second.y - another.y) > 0)) {
      if (Consts::Equals(first.y, second.y)) {
        continue;
      }
      Point temp = intersection(Line(first, second), Line(another, 0));
      if (temp.x >= another.x) {
        ++count;
      }
    }
  }
  return count % 2 == 1;
};

std::pair<Line, Line> Ellipse::directrices() const {
  Line fst(1, 0, -pow(a(), 2) / c());
  Line snd(1, 0, -pow(a(), 2) / c());
  return {fst, snd};
}

void Ellipse::rotate(const Point& center, double angle) {
  first.rotate(center, angle);
  second.rotate(center, angle);
}

void Ellipse::reflect(const Point& center) {
  first.reflect(center);
  second.reflect(center);
}

void Ellipse::reflect(const Line& axis) {
  first.reflect(axis);
  second.reflect(axis);
}

void Ellipse::scale(const Point& center, double coefficient) {
  first.scale(center, coefficient);
  second.scale(center, coefficient);
  rad *= coefficient;
}

double Ellipse::perimeter() const {
  return M_PI * (3 *(a() + b()) - sqrt((3 * a() + b()) * (a() + 3 * b())));
}

double Ellipse::area() const { return M_PI *a() *b(); }

bool Ellipse::operator==(const Shape& another) const {
  const Ellipse* temp = dynamic_cast<const Ellipse* >(&another);
  return temp != nullptr &&
         ((first == temp->first && second == temp->second) ||
          (first == temp->second && second == temp->first))&& 
         Consts::Equals(rad, temp->rad);
}

bool Ellipse::operator!=(const Shape& another) const {
  return !(*this == another);
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  const Ellipse* temp = dynamic_cast<const Ellipse* >(&another);
  return temp != nullptr &&
         Consts::Equals(Consts::Distance(first, second),
                        Consts::Distance(temp->first, temp->second)) &&
         Consts::Equals(rad, temp->rad);
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  const Ellipse* temp = dynamic_cast<const Ellipse* >(&another);
  return temp != nullptr &&
         Consts::Equals(Consts::Distance(first, second) *temp->rad,
                        Consts::Distance(temp->first, temp->second)* 
                            temp->rad);
}

bool Ellipse::containsPoint(const Point& another) const {
  double dist =
      Consts::Distance(another, first) + Consts::Distance(another, second);
  return (dist < rad + Consts::kNull);
}

Rectangle::Rectangle(const Point& first, const Point& third,
                     double coefficient) {
  if (coefficient < 1) {
    coefficient = 1 / coefficient;
  }
  double diag = Consts::Distance(first, third);
  double small = sqrt(pow(diag, 2) / (1 + pow(coefficient, 2)));
  double fcord = (third.x - first.x) *small / diag + first.x;
  double scord = (third.y - first.y) *small / diag + first.y;
  Point second(fcord, scord);
  double angle = atan(coefficient) *Consts::kStraight_angle / M_PI;
  second.rotate(first, angle);
  Point center(Consts::Middle(first, third));
  Point fourth = second;
  fourth.reflect(center);
  vertices = {first, second, third, fourth};
}

std::pair<Line, Line> Rectangle::diagonals() const {
  Line diag1(vertices[0], vertices[2]);
  Line diag2(vertices[1], vertices[3]);
  return {diag1, diag2};
}

Circle Square::circumscribedCircle() const {
  Point O = center();
  double radius = side() / sqrt(2);
  return Circle(O, radius);
}

Circle Square::inscribedCircle() const {
  Point O = center();
  double radius = side() / 2;
  return Circle(O, radius);
}

Point Triangle::circumcenter() const {
  Point mid1 = Consts::Middle(vertices[0], vertices[1]);
  Line side1 = Line(vertices[0], vertices[1]);
  Line serp1 = Consts::Perpendicular(mid1, Line(side1));
  Point mid2 = Consts::Middle(vertices[1], vertices[2]);
  Line side2 = Line(vertices[1], vertices[2]);
  Line serp2 = Consts::Perpendicular(mid2, Line(side2));
  return intersection(serp1, serp2);
}

Circle Triangle::circumscribedCircle() const {
  return Circle(circumcenter(), Consts::Distance(circumcenter(), vertices[0]));
}

Line Triangle::bisector(int index) const {
  Point now = vertices[index];
  Point next = vertices[(index + 1) % 3];
  Point prev = vertices[(2 + index) % 3];
  double first = Consts::Distance(now, prev);
  double second = Consts::Distance(now, next);
  double xx = (prev.x *second + next.x *first) / (first + second);
  double yy = (prev.y *second + next.y *first) / (first + second);
  Point temp(xx, yy);
  return Line(now, temp);
}

Circle Triangle::inscribedCircle() const {
  Line bis1 = bisector(1);
  Line bis2 = bisector(2);
  Line side(vertices[0], vertices[1]);
  Point incenter = intersection(bis1, bis2);
  return Circle(incenter, Consts::Distance(incenter, side));
}

Point Triangle::centroid() const {
  Point first = Consts::Middle(vertices[0], vertices[1]);
  Point second = Consts::Middle(vertices[1], vertices[2]);
  return intersection(Line(first, vertices[2]), Line(second, vertices[0]));
}

Point Triangle::orthocenter() const {
  Line side1(vertices[0], vertices[1]);
  Line side2(vertices[1], vertices[2]);
  return intersection(Consts::Perpendicular(vertices[2], side1),
                      Consts::Perpendicular(vertices[0], side2));
}

Line Triangle::EulerLine() const { return Line(orthocenter(), centroid()); }

Circle Triangle::ninePointsCircle() const {
  Point center(Consts::Middle(circumcenter(), orthocenter()));
  return Circle(center, Consts::Distance(
                            center, Consts::Middle(vertices[0], vertices[1])));
}