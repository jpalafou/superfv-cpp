#pragma once
#include "array.hpp"
#include "mesh.hpp"
#include "boundary_conditions.hpp"

enum Integrator {forward_euler, ssprk2, ssprk3, rk4};

void update_fluxes(const Array& u, Array& f, Array& g, Array& h, int p, const BC_Params& bcparams, double gamma);

double compute_dt(const Array& u, double cfl, double gamma, const UniformMeshParams& meshparams);

void take_step(
    const Array& uold,
    Array& unew,
    Array& f,
    Array& g,
    Array& h,
    double dt,
    const UniformMeshParams& meshparams,
    const BC_Params& bcparams,
    int p,
    Integrator integrator,
    double gamma
);