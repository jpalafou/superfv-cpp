#include <iostream>
#include <vector>
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "common.hpp"
#include "hydro.hpp"
#include "initial_conditions.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "ode.hpp"

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
    SolverParams params{
        HydroParams{1.4, 0.8},
        UniformMeshParams{0.0, 1.0, 16},
        BC_Params{3, BC::periodic, BC::none},
        3,
        Integrator::ssprk3
    };
    Array u = superfv(Setup::SQUARE, params, 0.000001);

    std::cout << "u\n";
    print_1d_array(u);
}
