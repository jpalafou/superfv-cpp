#include <vector>
#include "array.hpp"

void stencil_sweep_1d(const Array& u_in, Array& u_out, const std::vector<float>& stencil) {
    size_t left_reach = stencil.size() / 2;
    size_t right_reach = (stencil.size() - 1) / 2;
    for (size_t v = 0; v < u_in.get_nvars(); v++) {
        for (size_t i = 0; i < u_in.get_nx(); i++) {
            // check that we aren't on a ghost cell
            if (i < left_reach || i >= u_in.get_nx() - right_reach) {continue;}

            // perform stencil multiplication
            double result = 0.0;
            for (size_t w = 0; w < stencil.size(); w++) {
                size_t i_t = i + w - left_reach;
                result += stencil[w] * u_in(v, i_t);
            }
            u_out(v, i) = result;
        }
    }
}

void stencil_sweep_2d(const Array& u_in, Array& u_out, const std::vector<float>& stencil, size_t axis) {
    if (axis != 1 && axis != 2) {throw std::invalid_argument("axis must be 1 or 2.");}
    size_t left_reach = stencil.size() / 2;
    size_t right_reach = (stencil.size() - 1) / 2;
    for (size_t v = 0; v < u_in.get_nvars(); v++) {
        for (size_t i = 0; i < u_in.get_nx(); i++) {
            for (size_t j = 0; j < u_in.get_ny(); j++) {
                // check that we aren't on a ghost cell
                if (axis == 1 && (i < left_reach || i >= u_in.get_nx() - right_reach)) {continue;}
                if (axis == 2 && (j < left_reach || j >= u_in.get_ny() - right_reach)) {continue;}

                // perform stencil multiplication
                double result = 0.0;
                for (size_t w = 0; w < stencil.size(); w++) {
                    size_t i_t = i, j_t = j;
                    if (axis == 1) {    i_t += w - left_reach;}
                    else {              j_t += w - left_reach;}
                    result += stencil[w] * u_in(v, i_t, j_t);
                }
                u_out(v, i, j) = result;
            }
        }
    }
}

void stencil_sweep_3d(const Array& u_in, Array& u_out, const std::vector<float>& stencil, size_t axis) {
    if (axis < 1 || axis > 3) {throw std::invalid_argument("axis must be 1, 2 or 3.");}
    size_t left_reach = stencil.size() / 2;
    size_t right_reach = (stencil.size() - 1) / 2;
    for (size_t v = 0; v < u_in.get_nvars(); v++) {
        for (size_t i = 0; i < u_in.get_nx(); i++) {
            for (size_t j = 0; j < u_in.get_ny(); j++) {
                for (size_t k = 0; k < u_in.get_nz(); k++) {
                    // check that we aren't on a ghost cell
                    if (axis == 1 && (i < left_reach || i >= u_in.get_nx() - right_reach)) {continue;}
                    if (axis == 2 && (j < left_reach || j >= u_in.get_ny() - right_reach)) {continue;}
                    if (axis == 3 && (k < left_reach || k >= u_in.get_nz() - right_reach)) {continue;}

                    // perform stencil multiplication
                    double result = 0.0;
                    for (size_t w = 0; w < stencil.size(); w++) {
                        size_t i_t = i, j_t = j, k_t = k;
                        if (axis == 1) {        i_t += w - left_reach;}
                        else if (axis == 2) {   j_t += w - left_reach;}
                        else {                  k_t += w - left_reach;}
                        result += stencil[w] * u_in(v, i_t, j_t, k_t);
                    }
                    u_out(v, i, j, k) = result;
                }
            }
        }
    }
}

void stencil_sweep(const Array& u_in, Array& u_out, const std::vector<float>& stencil, size_t axis) {
    if (!has_same_shape(u_in, u_out)) {throw std::invalid_argument("u_in and u_out must have same shape.");}
    if (u_in.get_ndim() == 2) {
        stencil_sweep_1d(u_in, u_out, stencil);
        return;
    }
    if (u_in.get_ndim() == 3) {
        stencil_sweep_2d(u_in, u_out, stencil, axis);
        return;
    }
    if (u_in.get_ndim() == 4) {
        stencil_sweep_3d(u_in, u_out, stencil, axis);
        return;
    }
}

