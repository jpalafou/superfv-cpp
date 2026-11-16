#pragma once
#include "array.hpp"
#include "stencils.hpp"

void stencil_sweep_1d(const Array& u_in, Array& u_out, const Stencil& stencil, size_t istart = 0);

void stencil_sweep(
    const Array& u_in,
    Array& u_out,
    const Stencil& stencil,
    size_t axis = 1,
    size_t istart = 0,
    size_t jstart = 0,
    size_t kstart = 0
);