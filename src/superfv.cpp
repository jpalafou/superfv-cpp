#include <iostream>
#include <vector>
#include "array.hpp"
#include "boundary_conditions.hpp"
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
    Array uold = ic(mesh, Setup::SQUARE);
    Array f(uold.get_nvars(), uold.get_nx() + 1);
    Array unew = get_empty_copy(uold);
    Array dummy = Array();

    double dt = 0.8 / (2.0 / mesh.get_hx());
    take_step(uold, unew, f, dummy, dummy, mesh.get_hx(), mesh.get_hy(), mesh.get_hz(), dt, 0, 2, Integrator::forward_euler);

    std::cout << "uold\n";
    print_1d_array(uold);
    std::cout << "f\n";
    print_1d_array(f);
    std::cout << "unew\n";
    print_1d_array(unew);
}
