#include "io.hpp"

#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "array.hpp"
#include "common.hpp"
#include "utils.hpp"

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

std::vector<Real> get_snapshot_schedule(Real tstop, std::size_t nout) {
    if (nout == 1) {return {tstop};}
    else if (nout > 1) {return linspace(0, tstop, nout);}
    return {};
}

void write_snapshot(const std::filesystem::path& output_dir, const Array& u, Real t, std::size_t n) {
    std::filesystem::create_directories(output_dir);
    const auto path = output_dir / std::format("snapshot_{}.bin", n);
    std::ofstream file(path, std::ios::binary);
    if (!file) {throw std::runtime_error("cannot open snapshot: " + path.string());}

    file.write(reinterpret_cast<const char*>(&t), sizeof(t));
    const auto size = u.get_nvars() * u.get_nx() * u.get_ny() * u.get_nz();
    if (size > 0) {
        file.write(
            reinterpret_cast<const char*>(&u(0, 0, 0, 0)),
            static_cast<std::streamsize>(size * sizeof(Real))
        );
    }
    file.close();
    if (!file) {throw std::runtime_error("cannot write snapshot: " + path.string());}
    std::cout << "> wrote a snapshot to " << path.string() << std::endl;
}
