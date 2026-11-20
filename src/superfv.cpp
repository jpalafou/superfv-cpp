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
    double gamma = 1.4;
    UniformMeshParams meshparams = UniformMeshParams{2, 0.0, 1.0, 8};
    UniformMesh mesh = UniformMesh(meshparams);
    Array uold = ic(mesh, Setup::SQUARE, gamma);
    Array f(uold.get_nvars(), uold.get_nx() + 1);
    Array unew = get_empty_copy(uold);
    Array dummy = Array();
    int p = 0;
    Integrator integrator = Integrator::forward_euler;

    double dt = compute_dt(uold, 0.8, gamma, meshparams);
    take_step(uold, unew, f, dummy, dummy, dt, meshparams, p, integrator, gamma);

    std::cout << "uold\n";
    print_1d_array(uold);
    std::cout << "f\n";
    print_1d_array(f);
    std::cout << "unew\n";
    print_1d_array(unew);
}
