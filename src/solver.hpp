#pragma once
#include <cstddef>
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "common.hpp"
#include "hydro.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"
#include "ode.hpp"

struct SolverParams{
    const HydroParams hydro;
    const UniformMeshParams mesh;
    const BC_Params bc;
    const std::size_t p;
    const Integrator integrator;
    const Setup setup;
    const Real tstop;
};

Array superfv(const SolverParams& params);
