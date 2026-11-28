#pragma once
#include <cstddef>
#include "array.hpp"
#include "boundary_conditions.hpp"
#include "common.hpp"
#include "hydro.hpp"
#include "initial_condition_defs.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "ode.hpp"
#include "schemes.hpp"

struct SolverParams{
    const HydroParams hydro;
    const UniformMeshParams mesh;
    const BC_Params bc;
    const IOParams io;
    const MOOD_Scheme scheme;
    const Integrator integrator;
    const Setup setup;
    const Real tstop;
    const std::size_t max_steps;
    const bool clip_last_dt;
};

Array superfv(const SolverParams& params);
