#include <iostream>
#include <vector>

#include "array.hpp"
#include "update.hpp"

int main() {
    const std::vector<double> values = {1.0, 0.8, 1.0, 0.2, 1.0, 0.5};
    Array u(1, values.size());
    Array f(1, values.size() - 1);

    for (std::size_t i = 0; i < values.size(); ++i) {
        u(0, i) = values[i];
    }

    compute_1d_fluxes(u, f, 1);

    std::cout << "Fluxes\n";
    for (std::size_t i = 0; i < f.get_nx(); ++i) {
        if (i > 0) {std::cout << ", ";}
        std::cout << f(0, i);
    }
    std::cout << '\n';
}
