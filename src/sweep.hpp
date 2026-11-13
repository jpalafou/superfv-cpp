#pragma once
#include "array.hpp"

void stencil_sweep(const Array& u_in, Array& u_out, const std::vector<float>& stencil, size_t axis = 1);