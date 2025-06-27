#include "../src/utilities.hpp"
#include "utilities.hpp"

#include <iostream>

template <typename T> void test(const std::string &name) {
  if constexpr (has_from_chars_v<T>) {
    std::cout << "from_chars(" << name << ") is supported\n";
  } else {
    std::cout << "from_chars(" << name << ") is NOT supported\n";
  }
}

int main() {
  test<bool>("bool");
  test<char>("char");
  test<int>("int");
  test<unsigned int>("unsigned int");
  test<long>("long");
  test<unsigned long>("unsigned long");
  test<long long>("long long");
  test<unsigned long long>("unsigned long long");
  test<float>("float");
  test<double>("double");
}