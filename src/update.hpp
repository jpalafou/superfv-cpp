#pragma once
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "common.hpp"
#include "mesh.hpp"

struct SolverParams;

enum Integrator {forward_euler, ssprk2, ssprk3, rk4};

void update_fluxes(const Array& u, Array& f, Array& g, Array& h, int p, const BC_Params& bcparams, Real gamma);

Real compute_dt(const Array& u, Real cfl, Real gamma, const UniformMeshParams& meshparams);

void take_step(const Array& uold, Array& unew, Array& f, Array& g, Array& h, Real dt, const SolverParams& params);