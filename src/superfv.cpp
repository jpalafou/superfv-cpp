#include <iostream>
#include <vector>

#include "array.hpp"
#include "sweep.hpp"

int main() {
    const std::vector<double> values = {1.0, 1.0, 2.0, 1.0, 1.0, 1.0};
    const Stencil stencil{{1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f}, 1};
    Array u_in(1, values.size());
    Array u_out(1, values.size());

    for (std::size_t i = 0; i < values.size(); ++i) {
        u_in(0, i) = values[i];
    }

    stencil_sweep(u_in, u_out, stencil);

    std::cout << "Input\tSweep\n";
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cout << u_in(0, i) << '\t' << u_out(0, i) << '\n';
    }
}
