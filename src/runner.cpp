#include <iostream>
#include "array.hpp"
#include "config.hpp"
#include "solver.hpp"

void print_1d_array(Array& arr) {
    std::cout <<         "rho: ";
    for (std::size_t i = 0; i < arr.get_nx(); ++i) {std::cout << arr(0, i) << ", ";}
    std::cout << "\n" << " vx: ";
    for (std::size_t i = 0; i < arr.get_nx(); ++i) {std::cout << arr(1, i) << ", ";}
    std::cout << "\n" << " vy: ";
    for (std::size_t i = 0; i < arr.get_nx(); ++i) {std::cout << arr(2, i) << ", ";}
    std::cout << "\n" << " vz: ";
    for (std::size_t i = 0; i < arr.get_nx(); ++i) {std::cout << arr(3, i) << ", ";}
    std::cout << "\n" << "  P: ";
    for (std::size_t i = 0; i < arr.get_nx(); ++i) {std::cout << arr(4, i) << ", ";}
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <problem.cfg>\n";
        return 1;
    }
    const SolverParams params = read_config(argv[1]);
    Array u = superfv(params);

    std::cout << "u\n";
    print_1d_array(u);
}
