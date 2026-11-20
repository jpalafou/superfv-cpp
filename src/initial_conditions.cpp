#include <stdexcept>
#include "array.hpp"
#include "hydro.hpp"
#include "initial_conditions.hpp"
#include "mesh.hpp"

Array get_array_from_mesh(const UniformMeshParams& meshparams) {
    const size_t nvars = 5;
    if      (meshparams.ndim() == 1) {return Array(nvars, meshparams.nx);}
    else if (meshparams.ndim() == 2) {return Array(nvars, meshparams.nx, meshparams.ny);}
    else if (meshparams.ndim() == 3) {return Array(nvars, meshparams.nx, meshparams.ny, meshparams.nz);}
    throw std::invalid_argument("ndim must be 1, 2, or 3.");
}

Array ic(const UniformMesh& mesh, Setup setup, double gamma, const std::vector<double>& params) {
    UniformMeshParams meshparams = mesh.get_params();
    Array out = get_array_from_mesh(meshparams);
    size_t ndim = meshparams.ndim();

    if (ndim == 1) {
        for (int i = 0; i < meshparams.nx; i++) {
            double rho = ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::rho, ndim, params, setup);
            double vx  = ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::vx,  ndim, params, setup);
            double vy  = ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::vy,  ndim, params, setup);
            double vz  = ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::vz,  ndim, params, setup);
            double P   = ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::P,   ndim, params, setup);
            Primitives    primitive{rho, vx, vy, vz, P};
            Conservatives conserved = prim_to_cons(primitive, gamma);
            out(0, i) = conserved.rho;
            out(1, i) = conserved.m1;
            out(2, i) = conserved.m2;
            out(3, i) = conserved.m3;
            out(4, i) = conserved.E;
        }
    } else {
        throw std::logic_error("2D or 3D ics not yet implemented.");
    }

    return out;
}