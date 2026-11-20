#pragma once
#include <vector>

struct UniformMeshParams {
    size_t    nghost;
    double     xlim0;
    double     xlim1;
    size_t        nx; 
    double     ylim0 = 0.0; 
    double     ylim1 = 1.0; 
    size_t        ny = 1;
    double     zlim0 = 0.0; 
    double     zlim1 = 1.0;
    size_t        nz = 1;

    double hx() const;
    double hy() const;
    double hz() const;
    size_t ndim() const;
};

class UniformMesh {
    public:
        UniformMesh(const UniformMeshParams params);
        const UniformMeshParams& get_params() const;
        const std::vector<double> x_centers() const;
        const std::vector<double> y_centers() const;
        const std::vector<double> z_centers() const;
        const std::vector<double> x_faces() const;
        const std::vector<double> y_faces() const;
        const std::vector<double> z_faces() const;
    private:
        UniformMeshParams params_;
        std::vector<double> x_centers_, y_centers_, z_centers_, x_faces_, y_faces_, z_faces_;
};