#pragma once
#include <vector>

class UniformMesh {
    public:
        UniformMesh(
            double xlim0,
            double xlim1,
            size_t    nx, 
            double ylim0 = 0.0, 
            double ylim1 = 1.0, 
            size_t    ny = 1, 
            double zlim0 = 0.0, 
            double zlim1 = 1.0,
            size_t    nz = 1
        );
        const std::vector<double> x_centers() const;
        const std::vector<double> y_centers() const;
        const std::vector<double> z_centers() const;
        const std::vector<double> x_faces() const;
        const std::vector<double> y_faces() const;
        const std::vector<double> z_faces() const;
        size_t get_nx() const;
        size_t get_ny() const;
        size_t get_nz() const;
        double get_hx() const;
        double get_hy() const;
        double get_hz() const;
        size_t get_ndim() const;
    private:
        double xlim0_, xlim1_, ylim0_, ylim1_, zlim0_, zlim1_;
        std::vector<double> x_centers_, y_centers_, z_centers_, x_faces_, y_faces_, z_faces_;
};