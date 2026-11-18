#pragma once
#include "array.hpp"
#include "boundary_conditions.hpp"

void compute_1d_fluxes(const Array& u, Array& f, int p);