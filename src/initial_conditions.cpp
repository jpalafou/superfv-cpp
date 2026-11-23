#include <iterator>
#include <stdexcept>
#include "array.hpp"
#include "common.hpp"
#include "hydro.hpp"
#include "initial_conditions.hpp"
#include "mesh.hpp"
#include "quadrature.hpp"

Array get_array_from_mesh(const UniformMeshParams& meshparams) {
    const size_t nvars = 5;
    if      (meshparams.ndim() == 1) {return Array(nvars, meshparams.nx);}
    else if (meshparams.ndim() == 2) {return Array(nvars, meshparams.nx, meshparams.ny);}
    else if (meshparams.ndim() == 3) {return Array(nvars, meshparams.nx, meshparams.ny, meshparams.nz);}
    throw std::invalid_argument("ndim must be 1, 2, or 3.");
}

Array ic_1d(const UniformMesh& mesh, size_t p, Setup setup, Real gamma, const std::vector<Real>& params) {
    UniformMeshParams meshparams = mesh.get_params();
    Real hx = meshparams.hx();
    Array out = get_array_from_mesh(meshparams);
    size_t ndim = meshparams.ndim();
    size_t ngl = p / 2 + 1;

    if (ngl == 1) {
        for (size_t i = 0; i < meshparams.nx; i++) {
            Primitives primitive_cell{
                ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::rho, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,  Variable::vx, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,  Variable::vy, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,  Variable::vz, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,   Variable::P, ndim, params, setup)
            };
            Conservatives conservative_cell = prim_to_cons(primitive_cell, gamma);
            out(0, i) = conservative_cell.rho;
            out(1, i) = conservative_cell.m1;
            out(2, i) = conservative_cell.m2;
            out(3, i) = conservative_cell.m3;
            out(4, i) = conservative_cell.E;
        }
        return out;
    }

    Quadrature quadrature = gauss_legendre::get(ngl);
    for (size_t i = 0; i < meshparams.nx; i++) {
        Real rho_ave = 0;
        Real  mx_ave = 0;
        Real  my_ave = 0;
        Real  mz_ave = 0;
        Real   E_ave = 0;
        for (size_t q = 0; q < ngl; q++) {
            Real x = mesh.x_centers()[i] + (quadrature.x[q] * hx / (Real{2}));
            Real w = quadrature.w[q] / Real{2};
            Primitives primitive_node{
                ic_func(mesh.x_centers()[i], 0.0, 0.0, Variable::rho, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,  Variable::vx, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,  Variable::vy, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,  Variable::vz, ndim, params, setup),
                ic_func(mesh.x_centers()[i], 0.0, 0.0,   Variable::P, ndim, params, setup)
            };
            Conservatives conservative_node = prim_to_cons(primitive_node, gamma);
            rho_ave += w * conservative_node.rho;
            mx_ave  += w * conservative_node.m1;
            my_ave  += w * conservative_node.m2;
            mz_ave  += w * conservative_node.m3;
            E_ave   += w * conservative_node.E;
        }
        out(0, i) = rho_ave;
        out(1, i) =  mx_ave;
        out(2, i) =  my_ave;
        out(3, i) =  mz_ave;
        out(4, i) =   E_ave;
    }

    return out;
}

Array ic(const UniformMesh& mesh, size_t p, Setup setup, Real gamma, const std::vector<Real>& params) {
    switch (mesh.get_params().ndim()) {
        case 1:
            return ic_1d(mesh, p, setup, gamma, params);
        case 2:
            throw std::invalid_argument("2D initial conditions not implemented yet.");
        case 3:
            throw std::invalid_argument("3D initial conditions not implemented yet.");
        default:
            throw std::invalid_argument("mesh must be 1D, 2D, or 3D.");
    }
}