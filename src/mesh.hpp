#pragma once
#include <cstddef>
#include <vector>
#include "common.hpp"

struct UniformMeshParams {
    const Real   xlim0;
    const Real   xlim1;
    const std::size_t nx;
    const Real   ylim0 = 0.0; 
    const Real   ylim1 = 1.0; 
    const std::size_t ny = 1;
    const Real   zlim0 = 0.0; 
    const Real   zlim1 = 1.0;
    const std::size_t nz = 1;

    Real hx() const;
    Real hy() const;
    Real hz() const;
    std::size_t ndim() const;
};

class UniformMesh {
    public:
        UniformMesh(const UniformMeshParams& params);
        const UniformMeshParams& get_params() const;
        const std::vector<Real>& x_centers() const;
        const std::vector<Real>& y_centers() const;
        const std::vector<Real>& z_centers() const;
        const std::vector<Real>& x_faces() const;
        const std::vector<Real>& y_faces() const;
        const std::vector<Real>& z_faces() const;
    private:
        UniformMeshParams params_;
        std::vector<Real> x_centers_, y_centers_, z_centers_, x_faces_, y_faces_, z_faces_;
};