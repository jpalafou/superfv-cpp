#include <stdexcept>
#include <algorithm>
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "constants.hpp"
#include "hydro.hpp"
#include "riemann_solver.hpp"
#include "stencils.hpp"
#include "sweep.hpp"
#include "update.hpp"

constexpr size_t nghost = 2;

void compute_1d_fluxes(const Array& upad, Array& f, int p) {
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

void compute_fluxes(const Array& upad, Array& f, int p) {
    if (upad.get_ndim()==2) {
        compute_1d_fluxes(upad, f, p);
    } else {
        throw std::logic_error("2D and 3D fluxes not implemented yet");
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

double compute_dt(const Array& u, double cfl, double gamma, double hx, double hy, double hz) {
    switch (u.get_ndim()) {
        case 2:
            return compute_1d_dt(u, cfl, gamma, hx);
        case 3:
            return compute_2d_dt(u, cfl, gamma, hx, hy);
        case 4:
            return compute_3d_dt(u, cfl, gamma, hx, hy, hz);
        default:
            throw std::invalid_argument("got invalid ndim.");
    }
}