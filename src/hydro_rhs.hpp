#pragma once
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "common.hpp"
#include "mesh.hpp"

void update_fluxes(const Array& u, Array& f, Array& g, Array& h, int p, const BC_Params& bcparams, Real gamma);

void update_dudt(Array& dudt, const Array& f, const Array& g, const Array& h, const UniformMeshParams& meshparams);

Real compute_dt(const Array& u, Real cfl, Real gamma, const UniformMeshParams& meshparams);