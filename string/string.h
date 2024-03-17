#include <iostream>
#include <algorithm>
#include <cstring>

class String {
private:
  size_t sz_ = 0;  
  size_t cap_ = 8;
  char* arr_ = nullptr;

  String(size_t count) : sz_(count), cap_(sz_ + 1), arr_(new char[cap_]) {}

  void change_cap(size_t len) {
    char* temp = new char[len + 1];
    memcpy(temp, arr_, sz_ + 1);
    std::swap(arr_, temp);
    delete[] temp;
    cap_ = len + 1;
  }

  size_t secret_find(const String& str, bool reverse) const {
    size_t begin = 0;
    size_t end = sz_ - str.sz_;
    int step = 1;

    if (reverse) {
      begin = sz_ - str.sz_;
      end = -1;
      step = -1;
    }

    for (size_t i = begin; i < end; i += step) { 
      bool good = true;

      for (size_t j = 0; j < str.sz_; ++j) {
        if (str[j] != arr_[i + j]) {
          good = false;
          break;
        }
      }
      if (good) {
        return i;
      }
    }

    return sz_;
  }

public:
  void swap(String& other) {
    std::swap(arr_, other.arr_);
    std::swap(sz_, other.sz_);
    std::swap(cap_, other.cap_);
  }

  String() : arr_(new char[cap_]) {
    arr_[sz_] = '\0';
  };

  String(const size_t& count,   char symbol) : String(count) {
    std::memset(arr_, symbol, count);
    arr_[sz_] = '\0';
  }

  String(const char* str) : String(strlen(str)) {
    std::memcpy(arr_, str, sz_);
    arr_[sz_] = '\0';
  }

  String(const char symbol) : String(1, symbol) {}

  String(const String& other): String(other.sz_) {
    memcpy(arr_, other.arr_, sz_);
    arr_[sz_] = '\0';
  }

  String& operator=(const String& other) {
    if (this == &other) {
      return *this;
    }
    delete[] arr_;
    arr_ = new char[other.cap_];
    sz_ = other.sz_;
    cap_ = other.cap_;
    memcpy(arr_, other.arr_, sz_ + 1);
    return *this;
  }

  ~String() { 
    delete[] arr_;
  }

  size_t size() const {
    return sz_;
  }

  size_t length() const {
    return sz_;
  }

  size_t capacity() const {
    return cap_ - 1;
  }

  bool empty() const {
    return sz_ == 0;
  }

  void clear() {
    sz_ = 0;
  }

  void shrink_to_fit() {
    change_cap(sz_);
  }

  void push_back(const char symbol) {
    if (sz_ + 1 == cap_) {
      change_cap(2 * sz_);
    }

    arr_[sz_] = symbol;
    sz_ += 1;
    arr_[sz_] = '\0';
  }

  void pop_back() {
    sz_ -= 1;
    arr_[sz_] = '\0';
  }

  char* data() {
    return arr_;
  }

  const char* data() const {
    return arr_;
  }

  char& operator[](size_t pos) {
    return arr_[pos];
  }

  const char& operator[](size_t pos) const {
    return arr_[pos];
  }

  char& front() {
    return arr_[0];
  }

  const char& front() const {
    return arr_[0];
  }

  char& back() {
    return arr_[sz_ - 1];
  }

  const char& back() const {
    return arr_[sz_ - 1];
  }

  size_t find(const String& str) const {
    return secret_find(str, false);
  }

  size_t rfind(const String& str) const {
    return secret_find(str, true);
  }

  String& operator+=(const String& other) {
    if (cap_ < sz_ + other.sz_ + 1) {
      change_cap(2 * (sz_ + other.sz_));
    }

    memcpy((*this).arr_ + sz_, other.arr_, other.sz_ + 1);
    sz_ += other.sz_;
    return *this;
  }

  String substr(size_t start, size_t count) const {
    size_t len = std::min(sz_ - start, count);

    if (sz_ < start) {
      len = 0;
    }

    String result(len);
    memcpy(result.arr_, arr_ + start, len);
    result[len] = '\n';
    return result;
  }
  friend bool operator==(const String&, const String&);
};

std::ostream& operator<<(std::ostream& out, const String& str) {
  out << str.data();
  return out;
}

std::istream& operator>>(std::istream& in, String& str) {
  str.clear();
  char symbol;
  in.get(symbol);

  while (!isspace(symbol) && !in.eof()) {
    str.push_back(symbol);
    in.get(symbol);
  }

  return in;
}

String operator+(const String& str1, const String& str2) {
  String result = str1;
  result += str2;
  return result;
}

bool operator==(const String& str1, const String& str2) {
  if (str1.size() != str2.size()) {
    return false;
  }

  for (size_t i = 0; i < str1.size(); ++i) {
    if (str1[i] != str2[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const String& str1, const String& str2) {
  return !(str1 == str2);
}

bool operator<(const String& str1, const String& str2) {
  for (size_t i = 0; i < std::min(str1.size(), str2.size()); ++i) {
    if (str1[i] > str2[i]) {
      return false;
    }
    else if (str1[i] < str2[i]) {
      return true;
    }
  }
  if (str1.size() < str2.size()) {
    return true;
  }
  return false;
}

bool operator>=(const String& str1, const String& str2) {
  return !(str1 < str2);
}

bool operator>(const String& str1, const String& str2) {
  return str2 < str1;
}

bool operator<=(const String& str1, const String& str2) {
  return !(str2 < str1);
}
