#include <stdexcept>
#include "array.hpp"
#include "stencils.hpp"

void stencil_sweep_1d(const Array& u_in, Array& u_out, const Stencil& stencil, size_t istart = 0) {
    const size_t istop = istart + u_out.get_nx();
    const size_t reach = stencil.left_reach();
    for (size_t v = 0; v < u_in.get_nvars(); v++) {
        for (size_t i = istart; i < istop; i++) {
            double result = 0.0;
            for (size_t w = 0; w < stencil.size(); w++) {
                size_t i_t = i + w - reach;
                result += stencil[w] * u_in(v, i_t);
            }
            u_out(v, i - istart) = result;
        }
    }
}

void stencil_sweep_2d(
    const Array& u_in,
    Array& u_out,
    const Stencil& stencil,
    size_t axis,
    size_t istart = 0,
    size_t jstart = 0
) {
    if (axis != 1 && axis != 2) {throw std::invalid_argument("axis must be 1 or 2.");}
    const size_t istop = istart + u_out.get_nx();
    const size_t jstop = jstart + u_out.get_ny();
    const size_t reach = stencil.left_reach();
    for (size_t v = 0; v < u_in.get_nvars(); v++) {
        for (size_t i = istart; i < istop; i++) {
            for (size_t j = jstart; j < jstop; j++) {
                double result = 0.0;
                for (size_t w = 0; w < stencil.size(); w++) {
                    size_t i_t = i, j_t = j;
                    if (axis == 1) {    i_t += w - reach;}
                    else {              j_t += w - reach;}
                    result += stencil[w] * u_in(v, i_t, j_t);
                }
                u_out(v, i - istart, j - jstart) = result;
            }
        }
    }
}

void stencil_sweep_3d(
    const Array& u_in,
    Array& u_out,
    const Stencil& stencil,
    size_t axis,
    size_t istart = 0,
    size_t jstart = 0,
    size_t kstart = 0
) {
    if (axis < 1 || axis > 3) {throw std::invalid_argument("axis must be 1, 2 or 3.");}
    const size_t istop = istart + u_out.get_nx();
    const size_t jstop = jstart + u_out.get_ny();
    const size_t kstop = kstart + u_out.get_nz();
    const size_t reach = stencil.left_reach();
    for (size_t v = 0; v < u_in.get_nvars(); v++) {
        for (size_t i = istart; i < istop; i++) {
            for (size_t j = jstart; j < jstop; j++) {
                for (size_t k = kstart; k < kstop; k++) {
                    double result = 0.0;
                    for (size_t w = 0; w < stencil.size(); w++) {
                        size_t i_t = i, j_t = j, k_t = k;
                        if (axis == 1) {        i_t += w - reach;}
                        else if (axis == 2) {   j_t += w - reach;}
                        else {                  k_t += w - reach;}
                        result += stencil[w] * u_in(v, i_t, j_t, k_t);
                    }
                    u_out(v, i - istart, j - jstart, k - kstart) = result;
                }
            }
        }
    }
}

void stencil_sweep(
    const Array& u_in,
    Array& u_out,
    const Stencil& stencil,
    size_t axis,
    size_t istart,
    size_t jstart,
    size_t kstart
) {
    if (!has_same_ndim(u_in, u_out)) {throw std::invalid_argument("u_in and u_out must have same ndim.");}
    if (u_in.get_nvars() != u_out.get_nvars()) {throw std::invalid_argument("u_in and u_out must have same nvars.");}
    if (u_in.get_ndim() == 2) {
        stencil_sweep_1d(u_in, u_out, stencil, istart);
        return;
    }
    if (u_in.get_ndim() == 3) {
        stencil_sweep_2d(u_in, u_out, stencil, axis, istart, jstart);
        return;
    }
    if (u_in.get_ndim() == 4) {
        stencil_sweep_3d(u_in, u_out, stencil, axis, istart, jstart, kstart);
        return;
    }
}

