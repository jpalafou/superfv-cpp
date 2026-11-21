#include <stdexcept>
#include <algorithm>
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "hydro.hpp"
#include "mesh.hpp"
#include "riemann_solver.hpp"
#include "stencils.hpp"
#include "sweep.hpp"
#include "update.hpp"

void update_1d_fluxes(const Array& upad, Array& f, int p, size_t nghost, double gamma) {
    Array ul = get_empty_copy(upad);
    Array ur = get_empty_copy(upad);
    Stencil  left = conservative_interpolation(p, Position::left);
    Stencil right = conservative_interpolation(p, Position::right);
    size_t reach = left.left_reach();

    // interpolate left and right faces
    for (size_t i = reach; i < upad.get_nx() - reach; i++) {
        Array interp1(upad.get_nvars(), 1), interp2(upad.get_nvars(), 1);
        stencil_sweep(upad, interp1, left, 1, i);
        stencil_sweep(upad, interp2, right, 1, i);
        for (size_t v = 0; v < upad.get_nvars(); v++) {
            ul(v, i) = interp1(v, 0);
            ur(v, i) = interp2(v, 0);
        }
    }

    //       0   1   2   3
    //       f   f   f   f
    // u   u | u   u   u   u
    // 0   1   2   3   4   5

    for (size_t i = 0; i < f.get_nx(); i++) {
        size_t iadj = i + nghost - 1;
        Conservatives interp1{ur(0, iadj), ur(1, iadj), ur(2, iadj), ur(3, iadj), ur(4, iadj)};
        iadj++;
        Conservatives interp2{ul(0, iadj), ul(1, iadj), ul(2, iadj), ul(3, iadj), ul(4, iadj)};
        Primitives wl = cons_to_prim(interp1, gamma);
        Primitives wr = cons_to_prim(interp2, gamma);
        Conservatives numflux = llf(wl, wr, gamma);

        f(0, i) = numflux.rho;
        f(1, i) = numflux.m1;
        f(2, i) = numflux.m2;
        f(3, i) = numflux.m3;
        f(4, i) = numflux.E;
    }
}

void update_fluxes(const Array& u, Array& f, Array& g, Array& h, int p, const BC_Params& bcparams, double gamma) {
    Array upad = get_padded_u(u, bcparams.nghost);
    apply_bc(u, upad, bcparams);

    switch (u.get_ndim()) {
        case 2:
            update_1d_fluxes(upad, f, p, bcparams.nghost, gamma);
            return;
        case 3:
            throw std::invalid_argument("2D fluxes not yet implemented.");
            return;
        case 4:
            throw std::invalid_argument("3D fluxes not yet implemented.");
            return;
        default:
            throw std::invalid_argument("upad must have ndim 2, 3, or 4.");
    }
}

double compute_1d_dt(const Array& u, double cfl, double gamma, double hx) {
    double max_speed_over_length = 0.0;
    for (size_t i = 0; i < u.get_nx(); i++) {
        Primitives w = Primitives{u(0, i), u(1, i), u(2, i), u(3, i), u(4, i)};
        double s = prim_to_cs(w, gamma);
        s += std::abs(w.v1);
        max_speed_over_length = std::max(s / hx, max_speed_over_length);
    }
    if (max_speed_over_length == 0.0) {throw std::runtime_error("cannot compute hydro dt with 0 velocity.");}
    double dt = cfl / max_speed_over_length;
    return dt;
}

double compute_2d_dt(const Array& u, double cfl, double gamma, double hx, double hy) {
    double max_speed_over_length_sum = 0.0;
    for (size_t i = 0; i < u.get_nx(); i++) {
        for (size_t j = 0; i < u.get_ny(); j++) {
            Primitives w = Primitives{u(0, i, j), u(1, i, j), u(2, i, j), u(3, i, j), u(4, i, j)};
            double cs = prim_to_cs(w, gamma);
            double speed_over_length_sum = (std::abs(w.v1) + cs) / hx;
            speed_over_length_sum       += (std::abs(w.v2) + cs) / hy;
            if (speed_over_length_sum > max_speed_over_length_sum) {
                max_speed_over_length_sum = speed_over_length_sum;
            }
        }
    }
    if (max_speed_over_length_sum == 0.0) {throw std::runtime_error("cannot compute hydro dt with 0 velocity.");}
    double dt = cfl / max_speed_over_length_sum;
    return dt;
}

double compute_3d_dt(const Array& u, double cfl, double gamma, double hx, double hy, double hz) {
    double max_speed_over_length_sum = 0.0;
    for (size_t i = 0; i < u.get_nx(); i++) {
        for (size_t j = 0; j < u.get_ny(); j++) {
            for (size_t k = 0; k < u.get_nz(); k++) {
                double rho = u(0, i, j, k);
                double  mx = u(1, i, j, k);
                double  my = u(2, i, j, k);
                double  mz = u(3, i, j, k);
                double   P = u(4, i, j, k);
                Primitives w = Primitives{rho, mx, my, mz, P};
                double cs = prim_to_cs(w, gamma);
                double speed_over_length_sum = (std::abs(w.v1) + cs) / hx;
                speed_over_length_sum       += (std::abs(w.v2) + cs) / hy;
                speed_over_length_sum       += (std::abs(w.v3) + cs) / hz;
                if (speed_over_length_sum > max_speed_over_length_sum) {
                    max_speed_over_length_sum = speed_over_length_sum;
                }
            }
        }
    }
    if (max_speed_over_length_sum == 0.0) {throw std::runtime_error("cannot compute hydro dt with 0 velocity.");}
    double dt = cfl / max_speed_over_length_sum;
    return dt;
}

double compute_dt(const Array& u, double cfl, double gamma, const UniformMeshParams& meshparams) {
    switch (u.get_ndim()) {
        case 2:
            return compute_1d_dt(u, cfl, gamma, meshparams.hx());
        case 3:
            return compute_2d_dt(u, cfl, gamma, meshparams.hx(), meshparams.hy());
        case 4:
            return compute_3d_dt(u, cfl, gamma, meshparams.hx(), meshparams.hy(), meshparams.hz());
        default:
            throw std::invalid_argument("got invalid ndim.");
    }
}

void update_dudt_1d(Array& dudt, const Array& f, double hx) {
    for (size_t v = 0; v < dudt.get_nvars(); v++) {
        for (size_t i = 0; i < dudt.get_nx(); i++) {
            dudt(v, i) = -(f(v, i+1) - f(v, i)) / hx;
        }
    }
    return;
}

void update_dudt_2d(Array& dudt, const Array& f, const Array& g, double hx, double hy) {
    for (size_t v = 0; v < dudt.get_nvars(); v++) {
        for (size_t i = 0; i < dudt.get_nx(); i++) {
            for (size_t j = 0; j < dudt.get_ny(); j++) {
                dudt(v, i, j)  = -(f(v, i+1, j) - f(v, i, j)) / hx;
                dudt(v, i, j) += -(g(v, i, j+1) - g(v, i, j)) / hy;
            }
        }
    }
    return;
}

void update_dudt_3d(Array& dudt, const Array& f, const Array& g, const Array& h, double hx, double hy, double hz) {
    for (size_t v = 0; v < dudt.get_nvars(); v++) {
        for (size_t i = 0; i < dudt.get_nx(); i++) {
            for (size_t j = 0; j < dudt.get_ny(); j++) {
                for (size_t k = 0; k < dudt.get_nz(); k++) {
                    dudt(v, i, j, k)  = -(f(v, i+1, j, k) - f(v, i, j, k)) / hx;
                    dudt(v, i, j, k) += -(g(v, i, j+1, k) - g(v, i, j, k)) / hy;
                    dudt(v, i, j, k) += -(h(v, i, j, k+1) - h(v, i, j, k)) / hz;
                }
            }
        }
    }
    return;
}

void update_dudt(Array& dudt, const Array& f, const Array& g, const Array& h, const UniformMeshParams& meshparams) {
    switch (dudt.get_ndim()) {
        case 2:
            update_dudt_1d(dudt, f, meshparams.hx());
            return;
        case 3:
            update_dudt_2d(dudt, f, g, meshparams.hx(), meshparams.hy());
            return;
        case 4:
            update_dudt_3d(dudt, f, g, h, meshparams.hx(), meshparams.hy(), meshparams.hz());
            return;
        default:
            throw std::invalid_argument("uold must have ndim 2, 3, or 4.");
    }
}

void update_unew(const Array& uold, const Array& dudt, Array& unew, double dt) {
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

void forward_euler_step(
    const Array& uold,
    Array& unew,
    Array& f, 
    Array& g, 
    Array& h, 
    double dt,
    const UniformMeshParams& meshparams,
    const BC_Params& bcparams,
    int p,
    double gamma
) {
    Array dudt = get_empty_copy(unew);
    update_fluxes(uold, f, g, h, p, bcparams, gamma);
    update_dudt(dudt, f, g, h, meshparams);
    update_unew(uold, dudt, unew, dt);
    return;
}

void ssprk2_step(
    const Array& uold,
    Array& unew,
    Array& f, 
    Array& g, 
    Array& h, 
    double dt,
    const UniformMeshParams& meshparams,
    const BC_Params& bcparams,
    int p,
    double gamma
) {
    Array dudt = get_empty_copy(unew);
    Array utmp = get_empty_copy(unew);

    update_fluxes(uold, f, g, h, p, bcparams, gamma);
    update_dudt(dudt, f, g, h, meshparams);
    update_unew(uold, dudt, utmp, dt);

    update_fluxes(utmp, f, g, h, p, bcparams, gamma);
    update_dudt(dudt, f, g, h, meshparams);
    
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

void ssprk3_step(
    const Array& uold,
    Array& unew,
    Array& f, 
    Array& g, 
    Array& h, 
    double dt,
    const UniformMeshParams& meshparams,
    const BC_Params& bcparams,
    int p,
    double gamma
) {
    Array dut0 = get_empty_copy(unew);
    Array dut1 = get_empty_copy(unew);
    Array dut2 = get_empty_copy(unew);
    Array utmp = get_empty_copy(unew);

    update_fluxes(uold, f, g, h, p, bcparams, gamma);
    update_dudt(dut0, f, g, h, meshparams);
    update_unew(uold, dut0, utmp, dt);

    update_fluxes(utmp, f, g, h, p, bcparams, gamma);
    update_dudt(dut1, f, g, h, meshparams);
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

    update_fluxes(utmp, f, g, h, p, bcparams, gamma);
    update_dudt(dut2, f, g, h, meshparams);

    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    double dudt_sum_term = dut0(v, i, j, k) + dut1(v, i, j, k) + (2.0 * dut2(v, i, j, k));
                    unew(v, i, j, k)  = uold(v, i, j, k) + (dt * dudt_sum_term / 6.0);
                }
            }
        }
    }
    return;
}

void rk4_step(
    const Array& uold,
    Array& unew,
    Array& f, 
    Array& g, 
    Array& h, 
    double dt,
    const UniformMeshParams& meshparams,
    const BC_Params& bcparams,
    int p,
    double gamma
) {
    Array dut0 = get_empty_copy(unew);
    Array dut1 = get_empty_copy(unew);
    Array dut2 = get_empty_copy(unew);
    Array dut3 = get_empty_copy(unew);
    Array utmp = get_empty_copy(unew);

    update_fluxes(uold, f, g, h, p, bcparams, gamma);
    update_dudt(dut0, f, g, h, meshparams);
    update_unew(uold, dut0, utmp, 0.5 * dt);

    update_fluxes(utmp, f, g, h, p, bcparams, gamma);
    update_dudt(dut1, f, g, h, meshparams);
    update_unew(uold, dut1, utmp, 0.5 * dt);

    update_fluxes(utmp, f, g, h, p, bcparams, gamma);
    update_dudt(dut2, f, g, h, meshparams);
    update_unew(uold, dut2, utmp, dt);

    update_fluxes(utmp, f, g, h, p, bcparams, gamma);
    update_dudt(dut3, f, g, h, meshparams);
    
    for (size_t v = 0; v < uold.get_nvars(); v++) {
        for (size_t i = 0; i < uold.get_nx(); i++) {
            for (size_t j = 0; j < uold.get_ny(); j++) {
                for (size_t k = 0; k < uold.get_nz(); k++) {
                    double dudt_sum_term = dut0(v, i, j, k) + (2.0 * dut1(v, i, j, k));
                    dudt_sum_term += (2.0 * dut2(v, i, j, k)) + dut3(v, i, j, k);
                    unew(v, i, j, k)  = uold(v, i, j, k) + (dt * dudt_sum_term / 6.0);
                }
            }
        }
    }
    return;
}

void take_step(
    const Array& uold,
    Array& unew,
    Array& f,
    Array& g,
    Array& h,
    double dt,
    const UniformMeshParams& meshparams,
    const BC_Params& bcparams,
    int p,
    Integrator integrator,
    double gamma
) {
    switch (integrator) {
        case Integrator::forward_euler:
            forward_euler_step(uold, unew, f, g, h, dt, meshparams, bcparams, p, gamma);
            return;
        case Integrator::ssprk2:
            ssprk2_step(uold, unew, f, g, h, dt, meshparams, bcparams, p, gamma);
            return;
        case Integrator::ssprk3:
            ssprk3_step(uold, unew, f, g, h, dt, meshparams, bcparams, p, gamma);
            return;
        case Integrator::rk4:
            rk4_step(uold, unew, f, g, h, dt, meshparams, bcparams, p, gamma);
            return;
        default:
            throw std::invalid_argument("invalid integrator.");
    }
}