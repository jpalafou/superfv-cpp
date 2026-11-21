#include "array.hpp"
#include "common.hpp"
#include "hydro.hpp"
#include "initial_conditions.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"
#include "solver.hpp"

Array hydro_solver(Setup setup, const SolverParams& params, Real tstop) {
    UniformMesh mesh = UniformMesh(params.mesh);
    Array uold = ic(mesh, setup, params.hydro.gamma);
    Array f(uold.get_nvars(), uold.get_nx() + 1, uold.get_ny(), uold.get_nz());
    Array g, h;
    if (params.mesh.ndim() >= 2) {
        g = Array(uold.get_nvars(), uold.get_nx(), uold.get_ny() + 1, uold.get_nz());
    } 
    if (params.mesh.ndim() == 3) {
        h = Array(uold.get_nvars(), uold.get_nx(), uold.get_ny(), uold.get_nz() + 1);
    }
    Array unew = get_empty_copy(uold);
    Real t = 0;

    while (true) {
        bool stop = false;
        Real dt = compute_dt(uold, 0.8, params.hydro.gamma, params.mesh);
        if (t + dt > tstop) {dt = tstop - t; stop=true;};
        take_step(uold, unew, f, g, h, dt, params);
        t += dt;
        if (t >= tstop || stop) {break;}
        copy_contents(unew, uold);
    }

    return unew;
}