#include "io.hpp"

#include <cstddef>
#include <format>
#include <iostream>
#include "array.hpp"
#include "common.hpp"

StepSummary get_step_summary(const Array& u, size_t step_number, Real t) {
    StepSummary out{.step_number=step_number, .t=t};
    for (size_t i = 0; i < u.get_nx(); i++) {
        for (size_t j = 0; j < u.get_ny(); j++) {
            for (size_t k = 0; k < u.get_nz(); k++) {
                // update step summary
                if (u(0, i, j, k) > out.rhomax){
                    out.rhomax = u(0, i, j, k);
                };
                if (u(0, i, j, k) < out.rhomin){
                    out.rhomin = u(0, i, j, k);
                };
                out.Etotal += u(4, i, j, k);
            }
        }
    }
    return out;
}

void printer(const StepSummary& step_summary, Real Eref) {
    std::cout << std::format(
        "step: {:>6}, t: {:>8.2e}, rho max: {:>9.2e}, rho min: {:>9.2e}, E cons: {:>9.2e}\n",
        step_summary.step_number, step_summary.t, step_summary.rhomax,
        step_summary.rhomin, step_summary.Etotal - Eref
    );
}
