#include <stdexcept>
#include "array.hpp"
#include "stencils.hpp"
#include "sweep.hpp"
#include "update.hpp"

void compute_1d_fluxes(const Array& u, Array& f, int p) {
    Array wl(u.get_nvars(), u.get_nx());
    Array wr(u.get_nvars(), u.get_nx());
    Stencil  left = conservative_interpolation(p, Position::left);
    Stencil right = conservative_interpolation(p, Position::right);
    size_t reach = left.left_reach();

    // interpolate left and right faces
    for (size_t i = reach; i < u.get_nx() - reach; i++) {
        Array interp1(u.get_nvars(), 1), interp2(u.get_nvars(), 1);
        stencil_sweep(u, interp1, left, 1, i);
        stencil_sweep(u, interp2, right, 1, i);
        for (size_t v = 0; v < u.get_nvars(); v++) {
            wl(v, i) = interp1(v, 0);
            wr(v, i) = interp2(v, 0);
        }
    }

    for (size_t v = 0; v < f.get_nvars(); v++) {
        for (size_t i = 0; i < f.get_nx(); i++) {
            f(v, i) = wl(v, i + reach);
        }
    }
}
