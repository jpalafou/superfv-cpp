#pragma once
#include <vector>

struct UniformMeshParams {
    const double     xlim0;
    const double     xlim1;
    const size_t        nx; 
    const double     ylim0 = 0.0; 
    const double     ylim1 = 1.0; 
    const size_t        ny = 1;
    const double     zlim0 = 0.0; 
    const double     zlim1 = 1.0;
    const size_t        nz = 1;

    double hx() const;
    double hy() const;
    double hz() const;
    size_t ndim() const;
};

class UniformMesh {
    public:
        UniformMesh(const UniformMeshParams& params);
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