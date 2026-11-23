#include <stdexcept>
#include "array.hpp"
#include "common.hpp"
#include "hydro_rhs.hpp"
#include "solver.hpp"

void update_unew(const Array& uold, const Array& dudt, Array& unew, Real dt) {
    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    unew(v, i, j, k)  = uold(v, i, j, k) + (dt * dudt(v, i, j, k));
                }
            }
        }
    }
    return;
}

void forward_euler_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params) {
    Array dudt = get_empty_copy(unew);
    update_fluxes(uold, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dudt, f, g, h, params.mesh);
    update_unew(uold, dudt, unew, dt);
    return;
}

void ssprk2_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params) {
    Array dudt = get_empty_copy(unew);
    Array utmp = get_empty_copy(unew);

    update_fluxes(uold, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dudt, f, g, h, params.mesh);
    update_unew(uold, dudt, utmp, dt);

    update_fluxes(utmp, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dudt, f, g, h, params.mesh);
    
    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    unew(v, i, j, k)  = 0.5 * uold(v, i, j, k);
                    unew(v, i, j, k) += 0.5 * (utmp(v, i, j, k) + (dt * dudt(v, i, j, k)));
                }
            }
        }
    }
    return;
}

void ssprk3_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params) {
    Array dut0 = get_empty_copy(unew);
    Array dut1 = get_empty_copy(unew);
    Array dut2 = get_empty_copy(unew);
    Array utmp = get_empty_copy(unew);

    update_fluxes(uold, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut0, f, g, h, params.mesh);
    update_unew(uold, dut0, utmp, dt);

    update_fluxes(utmp, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut1, f, g, h, params.mesh);
    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    utmp(v, i, j, k)  = uold(v, i, j, k);
                    utmp(v, i, j, k) += 0.25 * dt * dut0(v, i, j, k);
                    utmp(v, i, j, k) += 0.25 * dt * dut1(v, i, j, k);
                }
            }
        }
    }

    update_fluxes(utmp, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut2, f, g, h, params.mesh);

    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    Real dudt_sum_term = dut0(v, i, j, k) + dut1(v, i, j, k) + (Real{2} * dut2(v, i, j, k));
                    unew(v, i, j, k)  = uold(v, i, j, k) + (dt * dudt_sum_term / Real{6});
                }
            }
        }
    }
    return;
}

void rk4_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params) {
    Array dut0 = get_empty_copy(unew);
    Array dut1 = get_empty_copy(unew);
    Array dut2 = get_empty_copy(unew);
    Array dut3 = get_empty_copy(unew);
    Array utmp = get_empty_copy(unew);

    update_fluxes(uold, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut0, f, g, h, params.mesh);
    update_unew(uold, dut0, utmp, 0.5 * dt);

    update_fluxes(utmp, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut1, f, g, h, params.mesh);
    update_unew(uold, dut1, utmp, 0.5 * dt);

    update_fluxes(utmp, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut2, f, g, h, params.mesh);
    update_unew(uold, dut2, utmp, dt);

    update_fluxes(utmp, f, g, h, params.p, params.bc, params.hydro.gamma);
    update_dudt(dut3, f, g, h, params.mesh);
    
    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    Real dudt_sum_term = dut0(v, i, j, k) + (Real{2} * dut1(v, i, j, k));
                    dudt_sum_term += (Real{2} * dut2(v, i, j, k)) + dut3(v, i, j, k);
                    unew(v, i, j, k)  = uold(v, i, j, k) + (dt * dudt_sum_term / Real{6});
                }
            }
        }
    }
    return;
}

void take_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params) {
    switch (params.integrator) {
        case Integrator::forward_euler:
            forward_euler_step(uold, unew, f, g, h, dt, params);
            return;
        case Integrator::ssprk2:
            ssprk2_step(uold, unew, f, g, h, dt, params);
            return;
        case Integrator::ssprk3:
            ssprk3_step(uold, unew, f, g, h, dt, params);
            return;
        case Integrator::rk4:
            rk4_step(uold, unew, f, g, h, dt, params);
            return;
        default:
            throw std::invalid_argument("invalid integrator.");
    }
}