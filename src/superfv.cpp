#include <iostream>
#include <vector>
#include "array.hpp"
#include "initial_conditions.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"
#include "update.hpp"

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

int main() {
    UniformMesh mesh = UniformMesh(0.0, 1.0, 8);
    Array u = ic(mesh, Setup::SQUARE, {2.0, 1.0, 1.0, 0.0, 0.0, 1.0});
    Array f(u.get_nvars(), u.get_nx() - 1);

    compute_1d_fluxes(u, f, 1);

    std::cout << "u\n";
    print_1d_array(u);
    std::cout << "f\n";
    print_1d_array(f);
}
