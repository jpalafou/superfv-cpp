#include <stdexcept>
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
