#include <algorithm>
#include <cstring>
#include <string>
#include <vector>


class BigInteger {
private:
  std::vector<int> digits{};
  bool isNegative = false;
  static const int kBase = 10000;
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

  BigInteger (int number) {
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
  friend BigInteger operator""_bi(const char* str);
};

BigInteger operator""_bi(unsigned long long numb) { return BigInteger(numb); }
BigInteger operator""_bi(const char* str) { return BigInteger(str); }


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
      : numerator(number.Abs()), denominator(1), isNegative(!number.sign()){}

  Rational(const BigInteger& first, const BigInteger& second)
      : numerator(first.Abs()),
        denominator(second.Abs()),
        isNegative(first.sign() != second.sign()) {
    Reduce();
  }

  Rational(int number)
    : numerator(abs(number)), denominator(1), isNegative(number < 0){}

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
    }
    else if (isNegative) {
      numerator = other.numerator * denominator - numerator;
      isNegative = (numerator < 0);
    }
    else {
      numerator -= other.numerator * denominator;
      isNegative = (numerator < 0);
    }
    numerator = numerator.Abs();
    denominator *= other.denominator;
    Reduce();
    Check_null();
    return *this;
  }

  Rational& operator-=(const Rational& other) {
    return *this += (-other);
  }

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
      if (begin == 0 &&  result[0] == '-') {
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
};

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