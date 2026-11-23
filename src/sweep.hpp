#pragma once
#include <cstddef>
#include "array.hpp"
#include "stencils.hpp"

void stencil_sweep_1d(const Array& u_in, Array& u_out, const Stencil& stencil, std::size_t istart = 0);

void stencil_sweep(
    const Array& u_in,
    Array& u_out,
    const Stencil& stencil,
    std::size_t axis = 1,
    std::size_t istart = 0,
    std::size_t jstart = 0,
    std::size_t kstart = 0
);