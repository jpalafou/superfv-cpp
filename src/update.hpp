#pragma once
#include "array.hpp"
#include "boundary_conditions.hpp"

void compute_1d_fluxes(const Array& upad, Array& f, int p);

void compute_fluxes(const Array& upad, Array& f, int p);

double compute_dt(const Array& u, double cfl, double gamma, double hx, double hy, double hz);