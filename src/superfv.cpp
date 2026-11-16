#include <iostream>
#include <vector>

#include "array.hpp"
#include "hydro.hpp"
#include "update.hpp"

int main() {
    const std::vector<double> rho_values = {1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0};
    const std::vector<double>  vx_values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const std::vector<double>   P_values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    Array u(5, rho_values.size());
    Array f(5, rho_values.size() - 1);

    for (std::size_t i = 0; i < rho_values.size(); ++i) {
        u(0, i) = rho_values[i];
        u(1, i) =  vx_values[i];
        u(4, i) =   P_values[i];
    }

    compute_1d_fluxes(u, f, 1);

    std::cout << "Fluxes\n";
    for (std::size_t v = 0; v < f.get_nvars(); ++v) {
        for (std::size_t i = 0; i < f.get_nx(); ++i) {
            if (i > 0) {std::cout << ", ";}
            std::cout << f(v, i);
        }
        std::cout << '\n';
    }
}
