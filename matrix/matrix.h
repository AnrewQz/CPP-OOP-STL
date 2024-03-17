#include <array>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

class BigInteger {
private:
  std::vector<int> digits{};
  bool isNegative = false;
  static const int kBase = 100000;
  static const int step = 4;
  BigInteger(const std::string& other) {
    isNegative = false;
    int end = -1;
    if (!isdigit(other[0])) {
      isNegative = true;
      end = 0;
    }
    for (int i = other.length() - 1; i > end; i -= step) {
      std::string temp;
      for (int j = std::max(0, end + step - i); j < step; ++j) {
        temp += other[i - step + j + 1];
      }
      digits.push_back(std::stoi(temp));
    }
  }

  static int Bin_search(const BigInteger& divident, const BigInteger& divider) {
    int left = 0;
    int right = kBase;
    while (right - left > 1) {
      int mid = (left + right) / 2;
      if (!(divident < (mid * divider))) {
        left = mid;
      } else {
        right = mid;
      }
    }
    return left;
  }

  void Check_null() {
    if (*this == 0 || -*this == 0) {
      isNegative = false;
    }
  }

  void Remove_zeros() {
    for (int i = digits.size() - 1; i > 0; --i) {
      if (digits[i] != 0) {
        break;
      }
      digits.pop_back();
    }
  }

public:
  BigInteger() { digits = {0}; }

  BigInteger(int number) {
    if (number == 0) {
      digits.push_back(0);
      return;
    }
    int copy = number;
    if (copy < 0) {
      isNegative = true;
      copy *= -1;
    }
    while (copy) {
      digits.push_back(copy % kBase);
      copy /= kBase;
    }
  }

  BigInteger(const BigInteger& other) = default;
  BigInteger& operator=(const BigInteger& other) = default;
  ~BigInteger() = default;

  bool sign() const { return !isNegative; }

  BigInteger operator-() const {
    BigInteger result = *this;
    if (!(result == 0)) {
      result.isNegative = !isNegative;
    }
    return result;
  }

  BigInteger Abs() const {
    BigInteger result = *this;
    result.isNegative = false;
    return result;
  }

  BigInteger& operator++() {
    *this += 1;
    return *this;
  }

  BigInteger& operator--() {
    *this -= 1;
    return *this;
  }

  BigInteger operator++(int) {
    BigInteger result = *this;
    ++*this;
    return result;
  }

  BigInteger operator--(int) {
    BigInteger result = *this;
    --*this;
    return result;
  }

  BigInteger& operator+=(const BigInteger& other) {
    if (other == 0) {
      return *this;
    }
    if (isNegative == other.isNegative) {
      for (size_t i = digits.size(); i < other.digits.size(); ++i) {
        digits.push_back(0);
      }
      int carry = 0;
      for (size_t i = 0; i < digits.size(); ++i) {
        if (i < other.digits.size()) {
          digits[i] += other.digits[i];
        }
        digits[i] += carry;
        carry = digits[i] / kBase;
        digits[i] %= kBase;
      }
      if (carry > 0) {
        digits.push_back(1);
      }
      Check_null();
      return *this;
    }
    return *this -= (-other);
  }

  BigInteger& operator-=(const BigInteger& other) {
    if (other == 0) {
      return *this;
    }
    if (isNegative == other.isNegative) {
      if (!isNegative) {
        if (!(*this < other)) {
          int carry = 0;
          for (size_t i = 0; i < digits.size(); ++i) {
            if (i < other.digits.size()) {
              digits[i] -= other.digits[i];
            }
            digits[i] -= carry;
            if (digits[i] < 0) {
              digits[i] += kBase;
              carry = 1;
            } else {
              carry = 0;
            }
          }
        } else {
          *this = -(other - *this);
        }
      } else {
        *this = -other - (-*this);
      }
      Remove_zeros();
      Check_null();
      return *this;
    }
    return *this += (-other);
  }

  BigInteger& operator*=(const BigInteger& other) {
    isNegative = !(isNegative == other.isNegative);
    std::vector<int> temp(digits.size() + other.digits.size() + 1, 0);
    for (size_t i = 0; i < other.digits.size(); ++i) {
      int carry = 0;
      for (size_t j = 0; j < digits.size(); ++j) {
        temp[i + j] += other.digits[i] * digits[j] + carry;
        carry = temp[i + j] / kBase;
        temp[i + j] %= kBase;
      }
      temp[i + digits.size()] += carry;
    }
    digits = temp;
    Remove_zeros();
    Check_null();
    return *this;
  }

  BigInteger& operator/=(const BigInteger& other) {
    BigInteger temp = other.Abs();
    *this = this->Abs();
    if (*this < temp) {
      *this = 0;
      return *this;
    }
    BigInteger result = 0;
    long long ind = digits.size() - 1;
    BigInteger current = 0;
    while (ind != -1) {
      current *= kBase;
      current += digits[ind];
      int quotient = Bin_search(current, temp);
      result = kBase * result;
      result += quotient;
      ind -= 1;
      current -= quotient * temp;
    }
    bool sign = !(isNegative == other.isNegative);
    *this = result;
    isNegative = sign;
    return *this;
  }

  BigInteger& operator%=(const BigInteger& other) {
    *this -= other * (*this / other);
    return *this;
  }

  std::string toString() const {
    std::string result;
    for (size_t i = 0; i < digits.size() - 1; ++i) {
      std::string str = std::to_string(digits[i]);
      result.insert(0, str);
      result.insert(0, 4 - str.length(), '0');
    }
    result.insert(0, std::to_string(digits[digits.size() - 1]));
    if (isNegative) {
      result.insert(0, "-");
    }
    return result;
  }

  explicit operator bool() const { return !(*this == 0); }
  explicit operator int() const { return stoi(toString()); }

  friend std::istream& operator>>(std::istream& , BigInteger& );
  friend bool operator<(const BigInteger& , const BigInteger& );
  friend bool operator==(const BigInteger& , const BigInteger& );
  friend BigInteger operator-(const BigInteger& , const BigInteger& );
  friend BigInteger operator*(const BigInteger& , const BigInteger& );
  friend BigInteger operator/(const BigInteger& , const BigInteger& );
  friend BigInteger operator""_bi(const char *str);
};

BigInteger operator""_bi(unsigned long long numb) { return BigInteger(numb); }
BigInteger operator""_bi(const char *str) { return BigInteger(str); }

std::istream& operator>>(std::istream& in, BigInteger& number) {
  std::string str;
  in >> str;
  number = BigInteger(str);
  return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
  out << number.toString();
  return out;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result += second;
  return result;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result -= second;
  return result;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result *= second;
  return result;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result /= second;
  return result;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
  BigInteger result = first;
  result %= second;
  return result;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
  if (first.digits.size() != second.digits.size() ||
      first.isNegative != second.isNegative) {
    return false;
  }
  for (size_t i = 0; i < first.digits.size(); ++i) {
    if (first.digits[i] != second.digits[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
  return !(first == second);
}

bool operator<(const BigInteger& first, const BigInteger& second) {
  if (first.isNegative == second.isNegative) {
    if (first.digits.size() < second.digits.size()) {
      return !first.isNegative;
    } else if (first.digits.size() > second.digits.size()) {
      return first.isNegative;
    }
    for (int i = first.digits.size() - 1; i >= 0; --i) {
      if (first.digits[i] < second.digits[i]) {
        return !first.isNegative;
      } else if (first.digits[i] > second.digits[i]) {
        return first.isNegative;
      }
    }
    return false;
  }
  return first.isNegative;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
  return second < first;
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
  return !(first < second);
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
  return !(second < first);
}

BigInteger Gcd(const BigInteger& first, const BigInteger& second) {
  BigInteger first_abs = first.Abs();
  BigInteger second_abs = second.Abs();
  while (first_abs && second_abs) {
    if (first_abs > second_abs) {
      first_abs %= second_abs;
    } else {
      second_abs %= first_abs;
    }
  }
  return first_abs + second_abs;
}

class Rational {

private:
  BigInteger numerator = 0;
  BigInteger denominator = 1;
  bool isNegative = false;
  static const int kDigits = 324;

  void Reduce() {
    BigInteger gcd = Gcd(numerator, denominator);
    numerator /= gcd;
    denominator /= gcd;
  }

  void Check_null() {
    if (numerator == 0) {
      isNegative = false;
    }
  }

public:
  Rational() {}

  Rational(const BigInteger& number)
      : numerator(number.Abs()), denominator(1), isNegative(!number.sign()) {}

  Rational(const BigInteger& first, const BigInteger& second)
      : numerator(first.Abs()), denominator(second.Abs()),
        isNegative(first.sign() != second.sign()) {
    Reduce();
  }

  Rational(int number)
      : numerator(abs(number)), denominator(1), isNegative(number < 0) {}

  Rational(const Rational& other) = default;
  Rational& operator=(const Rational& other) = default;
  ~Rational() = default;

  Rational operator-() const {
    Rational result = *this;
    if (!(numerator == 0)) {
      result.isNegative = !isNegative;
    }
    return result;
  }

  Rational& operator+=(const Rational& other) {
    numerator *= other.denominator;
    if (isNegative == other.isNegative) {
      numerator += other.numerator * denominator;
    } else if (isNegative) {
      numerator = other.numerator * denominator - numerator;
      isNegative = (numerator < 0);
    } else {
      numerator -= other.numerator * denominator;
      isNegative = (numerator < 0);
    }
    numerator = numerator.Abs();
    denominator *= other.denominator;
    Reduce();
    Check_null();
    return *this;
  }

  Rational& operator-=(const Rational& other) { return *this += (-other); }

  Rational& operator*=(const Rational& other) {
    isNegative = (!(isNegative == other.isNegative));
    numerator *= other.numerator;
    denominator *= other.denominator;
    Reduce();
    Check_null();
    return *this;
  }

  Rational& operator/=(const Rational& other) {
    isNegative = (!(isNegative == other.isNegative));
    numerator *= other.denominator;
    denominator *= other.numerator;
    Reduce();
    Check_null();
    return *this;
  }

  std::string toString() const {
    std::string result;
    if (isNegative) {
      result += "-";
    }
    result += numerator.Abs().toString();
    if (denominator.Abs() != 1) {
      result += "/" + denominator.Abs().toString();
    }
    return result;
  }

  std::string asDecimal(size_t precision = 0) const {
    std::string result;
    BigInteger numerator_temp = numerator.Abs();
    BigInteger denominator_temp = denominator.Abs();
    if (isNegative) {
      result += "-";
    }
    BigInteger begin = numerator_temp / denominator_temp;
    result += begin.toString();
    if (precision == 0) {
      if (begin == 0 && result[0] == '-') {
        result.erase(result.begin());
      }
      return result;
    }
    result += ".";
    BigInteger temp = numerator_temp % denominator_temp;
    for (size_t i = 0; i < precision; ++i) {
      temp *= 10;
      result += (temp / denominator_temp).toString();
      temp %= denominator_temp;
    }
    return result;
  }

  explicit operator double() const { return std::stod(asDecimal(kDigits)); }

  friend bool operator==(const Rational& first, const Rational& second);
  friend bool operator<(const Rational& first, const Rational& second);
  friend std::istream& operator>>(std::istream& in, Rational& numb);
};

std::istream& operator>>(std::istream& in, Rational& numb) {
  BigInteger numer;
  in >> numer;
  numb.numerator = numer;
  numb.denominator = 1;
  return in;
}

Rational operator+(const Rational& first, const Rational& second) {
  Rational result = first;
  result += second;
  return result;
}

Rational operator-(const Rational& first, const Rational& second) {
  Rational result = first;
  result -= second;
  return result;
}

Rational operator*(const Rational& first, const Rational& second) {
  Rational result = first;
  result *= second;
  return result;
}

Rational operator/(const Rational& first, const Rational& second) {
  Rational result = first;
  result /= second;
  return result;
}

bool operator==(const Rational& first, const Rational& second) {
  return (first.isNegative == second.isNegative &&
          first.numerator == second.numerator &&
          first.denominator == second.denominator);
}

bool operator!=(const Rational& first, const Rational& second) {
  return !(first == second);
}

bool operator<(const Rational& first, const Rational& second) {
  if (first.isNegative == second.isNegative) {
    if (first.numerator * second.denominator <
        first.denominator * second.numerator) {
      return !first.isNegative;
    }
    if (first.numerator * second.denominator >
        first.denominator * second.numerator) {
      return first.isNegative;
    }
    return false;
  }
  return first.isNegative;
}

bool operator>(const Rational& first, const Rational& second) {
  return second < first;
}

bool operator>=(const Rational& first, const Rational& second) {
  return !(first < second);
}

bool operator<=(const Rational& first, const Rational& second) {
  return !(first > second);
}

constexpr int sqrt(int n, int i = 1){
    return n == i ? n : (i * i < n ? sqrt(n, i + 1) : i);
}

template <int N, int D> struct isPrimeHelper {
  static const bool result =
      (N % D == 0 ? false : isPrimeHelper<N, D - 1>::result);
};

template <int N> struct isPrimeHelper<N, 1> {
  static const bool result = true;
};

template <int N> struct isPrime {
  static const bool result = isPrimeHelper<N, sqrt(N)>::result;
};

template <> struct isPrime<1> { static const bool result = false; };

template <int N> const bool is_prime = isPrime<N>::result;

template <size_t N> class Residue {
private:
  size_t value = 0;
  static void gcdex(int a, int b, int& x, int& y) {
    if (a == 0) {
      x = 0;
      y = 1;
      return;
    }
    int x1 = 0;
    int y1 = 0;
    gcdex(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
  }

  Residue<N> oppositeNumber() const {
    // static_assert(isPrime<N>::result, "Trying to devide not in Field");
    int x = 0;
    int y = 0;
    gcdex(value, int(N), x, y);
    x = (x % int(N) + N);
    return Residue<N>(x);
  }

public:
  Residue<N>() = default;
  Residue<N>(const int number) : value(number % int(N)) {}
  Residue<N>(const Residue<N>& other) = default;
  ~Residue<N>() = default;
  explicit operator int() const { return value; }
  size_t Value() const { return value; }
  Residue<N>& operator=(const Residue<N>& other) = default;
  Residue<N>& operator+=(const Residue<N>& other);
  Residue<N>& operator-=(const Residue<N>& other);
  Residue<N>& operator*=(const Residue<N>& other);
  Residue<N>& operator%=(const Residue<N>& other);
  Residue<N>& operator/=(const Residue<N>& other);
};

template <size_t N>
Residue<N>& Residue<N>::operator+=(const Residue<N>& other) {
  // std::cerr << " " << this->value << " + " << other.value << " = ";
  value = (value + other.value) % int(N);
  // std::cerr << this->value << " (" << N << ")\n";
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator-=(const Residue<N>& other) {
  // std::cerr << " " << this->value << " - " << other.value << " = ";
  value = (N + value - other.value) %  int(N);
  // std::cerr << this->value << " (" << N << ")\n";
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator*=(const Residue<N>& other) {
  // std::cerr << " " << this->value << " * " << other.value << " = ";
  value = (value * other.value) % int(N);
  // std::cerr << this->value << " (" << N << ")\n";
  return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator%=(const Residue<N>& other) {
  value %= int(other.value);
  return *this;
}

template <size_t N>
Residue<N> operator+(const Residue<N> first, const Residue<N> second) {
  Residue<N> result(first);
  result += second;
  return result;
}

template <size_t N>
Residue<N> operator-(const Residue<N> first, const Residue<N> second) {
  Residue<N> result(first);
  result -= second;
  return result;
}

template <size_t N>
Residue<N> operator*(const Residue<N> first, const Residue<N> second) {
  Residue<N> result(first);
  result *= second;
  return result;
}

template <size_t N>
Residue<N> operator%(const Residue<N> first, const Residue<N> second) {
  Residue<N> result(first);
  result %= second;
  return result;
}

template <size_t N>
Residue<N>& Residue<N>::operator/=(const Residue<N>& other) {
  // std::cerr << " " << this->value << " / " << other.value << " = ";
  Residue temp(other.oppositeNumber());
  *this *= temp;
  // std::cerr << this->value << " (" << N << ")\n";
  return *this;
}

template <size_t N>
Residue<N> operator/(const Residue<N> first, const Residue<N> second) {
  Residue<N> result(first);
  result /= second;
  return result;
}

template <size_t N>
bool operator==(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() == right.Value();
}

template <size_t N>
bool operator!=(const Residue<N>& left, const Residue<N>& right) {
  return !(left == right);
}

template <size_t M, size_t N, typename Field = Rational> class Matrix {
private:
  std::vector<std::vector<Field>> arr =
      std::vector<std::vector<Field>>(M, std::vector<Field>(N, Field(0)));
  void multitoolFoo(Field& det, Matrix<M, N, Field>& reversed) const;

public:
  Matrix() = default;
  Matrix(const std::vector<std::vector<Field>>& data);
  Matrix(std::initializer_list<std::initializer_list<Field>> list);
  Matrix(const Matrix<M, N, Field>& other) = default;
  Matrix& operator=(const Matrix<M, N, Field>& other) = default;
  Matrix<M, N, Field>& operator*=(const Field& other);
  Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other);
  Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other);
  Matrix<M, N, Field>& operator*=(const Matrix<N, N, Field>& other);
  const std::vector<Field>& operator[](const size_t index) const;
  std::vector<Field>& operator[](const size_t index);
  Matrix<N, M, Field> transposed() const;
  Field trace() const;
  Field det() const;
  size_t rank() const;
  Matrix<M, N, Field> inverted() const;
  Matrix<M, N, Field>& invert();
  std::array<Field, N> getRow(size_t index);
  std::array<Field, M> getColumn(size_t index);
};

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix(const std::vector<std::vector<Field>>& data) {
  for (int i = 0; i < data.size(); ++i) {
    for (int j = 0; j < data.size(); ++j) {
      arr[i][j] = data[i][j];
    }
  }
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix(std::initializer_list<std::initializer_list<Field>> list) {
  size_t i = 0, j = 0;
  for (auto& row: list) {
    for (auto& obj: row) {
      arr[i][j] = Field(obj);
      ++j;
    }
    ++i;
    j = 0;
  }
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(const Field& other) {
  // std::cerr << "Matr * Field";
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      arr[i][j] *= other;
    }
  }
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field& first, const Matrix<M, N, Field>& second) {
  Matrix<M, N, Field> result(second);
  result *= first;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& 
Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field>& other) {
  // std::cerr << "Matr + Matr";
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      arr[i][j] += other.arr[i][j];
    }
  }
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> first,
                              const Matrix<M, N, Field> second) {
  Matrix<M, N, Field> result(first);
  result += second;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& 
Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field>& other) {
  // std::cerr << "Matr - Matr";
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      arr[i][j] -= other.arr[i][j];
    }
  }
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> first,
                              const Matrix<M, N, Field> second) {
  Matrix<M, N, Field> result(first);
  result -= second;
  return result;
}

template <size_t M, size_t N, typename Field>
const std::vector<Field>& 
Matrix<M, N, Field>::operator[](const size_t index) const {
  return arr[index];
}

template <size_t M, size_t N, typename Field>
std::vector<Field>& Matrix<M, N, Field>::operator[](const size_t index) {
  return arr[index];
}

template <size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& first,
                               const Matrix<N, K, Field>& second) {
  // std::cerr << "Matr * Matr";
  Matrix<M, K, Field> result;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < N; ++k) {
        result[i][j] += first[i][k] * second[k][j];
      }
    }
  }
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& 
Matrix<M, N, Field>::operator*=(const Matrix<N, N, Field>& other) {
  *this = *this * other;
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
  // std::cerr << "Matr transpose";
  Matrix<N, M, Field> result;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      result[j][i] = (*this)[i][j];
    }
  }
  return result;
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const {
  // std::cerr << "Matr trace";
  static_assert(M == N, "taking a trace from non-square matrix");
  Field result(0);
  for (size_t i = 0; i < N; ++i) {
    result += (*this)[i][i];
  }
  return result;
}

template <size_t N, typename Field> Matrix<N, N, Field> identityMatrix() {
  // std::cerr << "Matr identity";
  Matrix<N, N, Field> result;
  for (size_t i = 0; i < N; ++i) {
    result[i][i] = Field(1);
  }
}

template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::multitoolFoo(Field& det,
                                       Matrix<M, N, Field>& inverted) const {
  // std::cerr << "multitool";
  static_assert(M == N, "you can take det/invert only square matrix");
  Matrix<M, 2 * M, Field> extended;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < M; ++j) {
      extended[i][j] = (*this)[i][j];
    }
    extended[i][i + M] = Field(1);
  }
  size_t curLevel = 0;
  size_t swapCount = 0;
  for (size_t i = 0; i < N; ++i) {
    Field mainValue(extended[curLevel][i]);
    for (size_t j = curLevel; j < M; ++j) {
      if (extended[j][i] != Field(0)) {
        mainValue = extended[j][i];
        std::swap(extended[j], extended[curLevel]);
        if (j != curLevel) {
          ++swapCount;
        }
        ++curLevel;
        break;
      }
    }
    det *= mainValue;
    if (mainValue != Field(0)) {
      for (size_t j = i; j < 2 * M; ++j) {
        extended[curLevel - 1][j] /= mainValue;
      }
      for (size_t j = curLevel; j < M; ++j) {
        Field coef = extended[j][i];
        for (size_t k = i; k < 2 * M; ++k) {
          extended[j][k] -= extended[curLevel - 1][k] * coef;
        }
      }
    }
  }
  det = (swapCount % 2 == 0) ? det : Field(-1) * det;
  size_t curColumn = N;
  for (size_t i = M - 1; i > 0; --i) {
    Field mainValue(0);
    for (size_t j = 0; j < curColumn; ++j) {
      if (extended[i][j] != Field(0)) {
        mainValue = extended[i][j];
        curColumn = j;
        break;
      }
    }
    if (mainValue != Field(0)) {
      for (size_t k = 0; k < i; ++k) {
        Field coef = extended[k][curColumn] / extended[i][curColumn];
        for (size_t j = curColumn; j < 2 * M; ++j) {
          extended[k][j] -= extended[i][j] * coef;
        }
      }
    }
  }
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < M; ++j) {
      inverted[i][j] = extended[i][j + M];
    }
  }
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const {
  // std::cerr << "det";
  Field result(1);
  Matrix<M, N, Field> temp;
  multitoolFoo(result, temp);
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::invert() {
  // std::cerr << "Invert";
  Field temp(0);
  Matrix<M, N, Field> result;
  multitoolFoo(temp, result);
  *this = result;
  return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
  Matrix<M, N, Field> result(*this);
  result.invert();
  return result;
}

template <size_t M, size_t N, typename Field>
std::array<Field, N> Matrix<M, N, Field>::getRow(size_t index) {
  // std::cerr << "getrow";
  std::array<Field, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = (*this)[index][i];
  }
  return result;
}

template <size_t M, size_t N, typename Field>
std::array<Field, M> Matrix<M, N, Field>::getColumn(size_t index) {
  // std::cerr << "getcolumn";
  std::array<Field, M> result;
  for (size_t i = 0; i < M; ++i) {
    result[i] = (*this)[i][index];
  }
  return result;
}

template <size_t M, size_t N, typename Field>
bool operator==(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
  // std::cerr << "matr == matr";
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (first[i][j] != second[i][j]) {
        return false;
      }
    }
  }
  return true;
}


template <size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
  // std::cerr << "rank";
  Matrix<M, N, Field> copy = *this;
  size_t curLevel = 0;
  size_t result = M;
  for (size_t i = 0; i < N; ++i) {
    Field mainValue(copy[curLevel][i]);
    for (size_t j = curLevel; j < M; ++j) {
      if (copy[j][i] != Field(0)) {
        mainValue = copy[j][i];
        std::swap(copy[j], copy[curLevel]);
        ++curLevel;
        break;
      }
    }
    if (mainValue != Field(0)) {
      for (size_t j = i; j < N; ++j) {
        copy[curLevel - 1][j] /= mainValue;
      }
      for (size_t j = curLevel; j < M; ++j) {
        Field coef = copy[j][i];
        for (size_t k = i; k < N; ++k) {
          copy[j][k] -= copy[curLevel - 1][k] * coef;
        }
      }
    }
  }
  for (size_t i = 0; i < M; ++i) {
    bool nullRow = true;
    for (size_t j = 0; j < N; ++j) {
      if (copy[i][j] != Field(0)) {
        nullRow = false;
        break;
      }
    }
    if (nullRow) {
      --result;
    }
  }
  return result;
}