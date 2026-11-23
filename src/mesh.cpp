#include "mesh.hpp"

#include <stdexcept>
#include <vector>
#include "common.hpp"

Real UniformMeshParams::hx()   const {return (xlim1 - xlim0) / static_cast<Real>(nx);}
Real UniformMeshParams::hy()   const {return (ylim1 - ylim0) / static_cast<Real>(ny);}
Real UniformMeshParams::hz()   const {return (zlim1 - zlim0) / static_cast<Real>(nz);}
size_t UniformMeshParams::ndim() const {if (nz > 1) {return 3;} else if (ny > 1) {return 2;} else {return 1;}}

std::vector<Real> linspace(Real lower, Real upper, size_t n){
    std::vector<Real> out;
    if (n < 2) {throw std::invalid_argument("n must be greater than 1.");}
    Real h = (upper - lower) / (n - 1);
    for(int i = 0; i < n - 1; ++i) {
        out.push_back(lower + h * i);
    }
    out.push_back(upper); 
    return out;
}

UniformMesh::UniformMesh(const UniformMeshParams& params) : params_(params) {
    x_faces_ = linspace(params.xlim0, params.xlim1, params.nx + 1);
    y_faces_ = linspace(params.ylim0, params.ylim1, params.ny + 1);
    z_faces_ = linspace(params.zlim0, params.zlim1, params.nz + 1);
    for (size_t i = 0; i < params.nx; i++) {x_centers_.push_back((x_faces_[i] + x_faces_[i+1]) / Real{2});}
    for (size_t j = 0; j < params.ny; j++) {y_centers_.push_back((y_faces_[j] + y_faces_[j+1]) / Real{2});}
    for (size_t k = 0; k < params.nz; k++) {z_centers_.push_back((z_faces_[k] + z_faces_[k+1]) / Real{2});}
}
const UniformMeshParams& UniformMesh::get_params() const {return params_;}
const std::vector<Real>& UniformMesh::x_centers() const {return x_centers_;}
const std::vector<Real>& UniformMesh::y_centers() const {return y_centers_;}
const std::vector<Real>& UniformMesh::z_centers() const {return z_centers_;}
const std::vector<Real>&   UniformMesh::x_faces() const {return   x_faces_;}
const std::vector<Real>&   UniformMesh::y_faces() const {return   y_faces_;}
const std::vector<Real>&   UniformMesh::z_faces() const {return   z_faces_;}