#pragma once
#include "array.hpp"
#include "common.hpp"

struct SolverParams;

enum Integrator {forward_euler, ssprk2, ssprk3, rk4};

void take_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params);