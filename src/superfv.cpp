#include <iostream>
#include "array.hpp"

int main() {
    Array array(1, 5);

    array(0, 0) = 1.0;

    for (std::size_t i = 0; i < 5; ++i) {std::cout << array(0, i) << ", ";}
    std::cout << std::endl;
}