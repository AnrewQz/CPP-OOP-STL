#include <iostream>
#include <iterator>

template <typename T>
class Deque {
private:

  const static int inner_cap_ = 32;
  int outer_cap_ = 0;
  size_t size_ = 0;
  T** arr_ = nullptr;
  int front_out_ = 0;
  int front_inn_ = 0;
  int back_out_ = 0;
  int back_inn_ = 0;

  void swap_(Deque<T>& other) {
    std::swap(outer_cap_, other.outer_cap_);
    std::swap(size_, other.size_);
    std::swap(arr_, other.arr_);
    std::swap(front_out_, other.front_out_);
    std::swap(front_inn_, other.front_inn_);
    std::swap(back_out_, other.back_out_);
    std::swap(back_inn_, other.back_inn_);
  }

  void resize() {
    T** new_arr = new T*[outer_cap_ * 3 + 2]();
    front_out_ += outer_cap_ + 1;
    back_out_ += outer_cap_ + 1;
    for (int i = 0; i < outer_cap_; ++i) {
      new_arr[outer_cap_ + 1 + i] = arr_[i];
    }
    for (int i = 0; i < 3 * outer_cap_ + 2; ++i) {
      if (new_arr[i] == nullptr) {
        new_arr[i] = reinterpret_cast<T*>(new char*[sizeof(T) * inner_cap_]);
      } 
    }
    delete[] arr_;
    outer_cap_ = outer_cap_ * 3 + 2;
    arr_ = new_arr;
  }

public:

  template <bool IsConst>
  class CommonIterator;

  using iterator = CommonIterator<false>;
  using const_iterator = CommonIterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() { return iterator(arr_ + front_out_, front_inn_); }
  iterator end() { return iterator(arr_ + back_out_, back_inn_); }

  const_iterator begin() const { return const_iterator(const_cast<const T**>(arr_ + front_out_), front_inn_); }
  const_iterator end() const { return const_iterator(const_cast<const T**>(arr_ + back_out_), back_inn_); }
  const_iterator cbegin() const { return const_iterator(const_cast<const T**>(arr_ + front_out_), front_inn_); }
  const_iterator cend() const { return const_iterator(const_cast<const T**>(arr_ + back_out_), back_inn_); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
  const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

  Deque() = default;
  Deque(int numb);
  Deque(const Deque<T>& other);
  Deque(int numb, const T& elem);
  Deque& operator=(const Deque<T>& other);
  ~Deque();
  size_t size() const;
  T& operator[](size_t index);
  const T& operator[](size_t index) const;
  T& at(size_t index);
  const T& at(size_t index) const;
  void push_back(const T& elem);
  void push_front(const T& elem);
  void pop_back();
  void pop_front();
  void insert(iterator it, const T& elem);
  void erase(iterator it);
};


template <typename T>
template <bool IsConst>
class Deque<T>::CommonIterator {
public:

  using difference_type = ptrdiff_t;
  using value_type = std::conditional_t<IsConst, const T, T>;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::random_access_iterator_tag;

private:
 pointer* outer_ = nullptr;
 int inner_ = 0;

public:

  CommonIterator() = delete;
  CommonIterator(pointer* ptr, int ind) : outer_(ptr), inner_(ind) {}
  CommonIterator(const CommonIterator& other) = default;
  CommonIterator& operator=(const CommonIterator& other) = default;
  ~CommonIterator() = default;

  CommonIterator& operator+=(ptrdiff_t x) {
    if (x >= 0) {
      outer_ += (x + inner_) / inner_cap_;
      inner_ = (x + inner_) % inner_cap_;
      return *this;
    }
    return *this -= -x;
  }

  CommonIterator& operator-=(ptrdiff_t x) {
    if (x >= 0) {
      outer_ -= ((x + inner_cap_) - inner_ - 1) / inner_cap_;
      inner_ = ((inner_ + inner_cap_) - (x % inner_cap_)) % inner_cap_;
      return *this;
    }
    return *this += -x;
  }

  CommonIterator& operator++() {
    return *this += 1;
  }

  CommonIterator& operator--() {
    return *this -= 1;
  }

  CommonIterator operator++(int) {
    auto result(*this);
    *this += 1;
    return result;
  }

  CommonIterator operator--(int) {
    auto result(*this);
    *this -= 1;
    return result;
  }

  CommonIterator operator+(ptrdiff_t x) const {
    auto result(*this);
    result += x;
    return result;
  }

  CommonIterator operator-(ptrdiff_t x) const {
    auto result(*this);
    result -= x;
    return result;
  }

  auto operator<=>(const CommonIterator &) const = default;


  ptrdiff_t operator-(const CommonIterator &other) const {
    return (inner_cap_ * (outer_ - other.outer_) + inner_) - other.inner_;
  }

  reference operator*() const {
    return *(*outer_ + inner_);
  }

  pointer operator->() const {
    return *outer_ + inner_;
  }

  reference operator[](ptrdiff_t index) const { return *(*this + index); }

  operator const_iterator() {
    return const_iterator(*this);
  }
};


template <typename T>
Deque<T>::Deque(int numb) : outer_cap_(numb / inner_cap_ + 1), size_(numb),
front_out_(0), front_inn_(0), back_out_(0), back_inn_(0) {
  int count = 0;
  try {
    arr_ = new T*[outer_cap_];
    for (int i = 0; i < outer_cap_; ++i) {
      arr_[i] = reinterpret_cast<T*>(new char*[sizeof(T) * inner_cap_]);
    }
    while (count != numb) {
      int i = count / inner_cap_;
      int j = count % inner_cap_;
      new(arr_[i] + j) T();
      ++count;
    }
    back_out_ = size_ / inner_cap_;
    back_inn_ = size_ % inner_cap_;
  } catch (...) {
    count--;
    while (count != -1) {
      int i = count / inner_cap_;
      int j = count % inner_cap_;
      (arr_[i] + j)->~T();
      --count;
    }
    for (int i = 0; i < numb / inner_cap_ + 1; ++i) {
      delete[] reinterpret_cast<char*>(arr_[i]);
    }
    delete[] arr_;
    throw;
  }
}

template <typename T>
Deque<T>::Deque(const Deque<T>& other) : outer_cap_(other.size_ / inner_cap_ + 1),
         size_(other.size_), front_out_(0), front_inn_(0), back_out_(0), back_inn_(0) {
  int count = 0;
  try {
    arr_ = new T*[outer_cap_];
    for (int i = 0; i < outer_cap_; ++i) {
      arr_[i] = reinterpret_cast<T*>(new char*[sizeof(T) * inner_cap_]);
    }
    auto it = other.cbegin();
    while (count != int(other.size_)) {
      int i = count / inner_cap_;
      int j = count % inner_cap_;
      new(arr_[i] + j) T(*it);
      ++it;
      ++count;
    }
    back_out_ = size_ / inner_cap_;
    back_inn_ = size_ % inner_cap_;
  } catch(...) {
    count--;
    while (count != -1) {
      int i = count / inner_cap_;
      int j = count % inner_cap_;
      (arr_[i] + j)->~T();
      --count;
    }
    for (int i = 0; i < outer_cap_; ++i) {
      delete[] reinterpret_cast<char*>(arr_[i]);
    }
    delete[] arr_;
    throw;
  }
}

template <typename T>
Deque<T>::Deque(int numb, const T& elem) : outer_cap_(numb / inner_cap_ + 1),
         size_(numb), front_out_(0), front_inn_(0), back_out_(0), back_inn_(0) {
  int count = 0;
  try {
    arr_ = new T*[outer_cap_];
    for (int i = 0; i < outer_cap_; ++i) {
      arr_[i] = reinterpret_cast<T*>(new char*[sizeof(T) * inner_cap_]);
    }
    while (count != numb) {
      int i = count / inner_cap_;
      int j = count % inner_cap_;
      new(arr_[i] + j) T(elem);
      ++count;
    }
    back_out_ = size_ / inner_cap_;
    back_inn_ = size_ % inner_cap_;
  } catch(...) {
    count--;
    while (count != -1) {
      int i = count / inner_cap_;
      int j = count % inner_cap_;
      (arr_[i] + j)->~T();
      --count;
    }
    for (int i = 0; i < outer_cap_; ++i) {
      delete[] reinterpret_cast<char*>(arr_[i]);
    }
    delete[] arr_;
    throw;
  }
}

template <typename T>
Deque<T>& Deque<T>::operator=(const Deque<T>& other) {
  if (this == &other) {
    return *this;
  }
  try {
    Deque<T> temp(other);
    swap_(temp);
    return *this;
  } catch(...) {
    return *this;
    throw;
  }
}

template <typename T>
Deque<T>::~Deque() {
  iterator it = begin();
  while (it != end()) {
    it->~T();
    ++it;
  }
  for (int i = 0; i < outer_cap_; ++i) {
    delete[] reinterpret_cast<char*>(arr_[i]);
  }
  delete[] arr_;
}

template <typename T>
size_t Deque<T>::size() const {
  return size_;
}

template <typename T>
T& Deque<T>::operator[](size_t index) {
  return *(begin() + index);
}

template <typename T>
const T& Deque<T>::operator[](size_t index) const {
  return *(cbegin() + index);
}

template <typename T>
T& Deque<T>::at(size_t index) {
  if (index >= size_ || index < 0) {
    throw std::out_of_range("index out of range");
  }
  return (*this)[index];
}

template <typename T>
const T& Deque<T>::at(size_t index) const {
  if (index >= size_ || index < 0) {
    throw std::out_of_range("index out of range");
  }
  return (*this)[index];
}

template <typename T>
void Deque<T>::push_back(const T& elem) {
  try {
    if (back_out_ == outer_cap_ && back_inn_ == 0) {
      resize();
    }
    new (arr_[back_out_] + back_inn_) T(elem);
    back_out_ += (++back_inn_) / inner_cap_;
    back_inn_ %= inner_cap_;
    ++size_;
  } catch(...) {
    throw;
  }
}

template <typename T>
void Deque<T>::push_front(const T& elem) {
  if (front_out_ == 0 && front_inn_ == 0) {
    resize();
  }
  ++size_;
  try {
    if (front_inn_ > 0) {
      new (arr_[front_out_] + --front_inn_) T(elem);
    } else {
      front_inn_ = inner_cap_ - 1;
      new (arr_[--front_out_] + front_inn_) T(elem);
    }
  } catch(...){
    if (front_inn_ != inner_cap_ - 1) {
      front_inn_++;
    } else {
      front_inn_ = 0;
      front_out_++;
    }
    throw;
  }
}

template <typename T>
void Deque<T>::pop_back() {
  back_out_ -= back_inn_ == 0 ? 1 : 0;
  back_inn_ = (back_inn_ - 1 + inner_cap_) % inner_cap_;
  arr_[back_out_][back_inn_].~T();
  --size_;
}

template <typename T>
void Deque<T>::pop_front() {
  arr_[front_out_][front_inn_].~T();
  front_out_ += front_inn_ == inner_cap_ - 1 ? 1 : 0;
  front_inn_ = (front_inn_ + 1) % inner_cap_;
  --size_;
}

template <typename T>
void Deque<T>::insert(iterator it, const T& elem) {
  T temp(elem);
  while (it != end()) {
    std::swap(*(it++), temp);
  }
  push_back(temp);
}

template <typename T>
void Deque<T>::erase(iterator it) {
  while (it != end()- 1) {
    std::swap(*it, *(++it));
  }
  pop_back();
}