#include <iostream>

#include "smart_pointers.h"
// #include "example2.h"
// #include "example.h"


struct Base {
  int x = 0;
};

struct Derived : Base {
  int y = 0;
};

int main() {
  int x = 1;
  SharedPtr<int>(&x, std::default_delete<int>(), std::allocator<int>());
}