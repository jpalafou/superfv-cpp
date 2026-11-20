#pragma once
#include "array.hpp"
#include "boundary_conditions.hpp"

enum Integrator {forward_euler, ssprk2, ssprk3, rk4};

void update_fluxes(const Array& u, Array& f, Array& g, Array& h, int p, size_t nghost);

double compute_dt(const Array& u, double cfl, double gamma, double hx, double hy, double hz);

void take_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, double hx, double hy, double hz, double dt, int p, size_t nghost, Integrator integrator);