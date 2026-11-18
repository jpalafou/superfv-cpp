#include <vector>
#include "mesh.hpp"

std::vector<double> linspace(double lower, double upper, size_t n){
    std::vector<double> out;
    if (n < 2) {throw std::invalid_argument("n must be greater than 1.");}
    double h = (upper - lower) / (n - 1);
    for(int i = 0; i < n - 1; ++i) {
        out.push_back(lower + h * i);
    }
    out.push_back(upper); 
    return out;
}

UniformMesh::UniformMesh(double xlim0, double xlim1, size_t nx, double ylim0, double ylim1, size_t ny, double zlim0, double zlim1, size_t nz) {
    xlim0_ = xlim0;
    xlim1_ = xlim1;
    ylim0_ = ylim0;
    ylim1_ = ylim1;
    zlim0_ = zlim0;
    zlim1_ = zlim1;
    x_faces_ = linspace(xlim0, xlim1, nx + 1);
    y_faces_ = linspace(ylim0, ylim1, ny + 1);
    z_faces_ = linspace(zlim0, zlim1, nz + 1);
    for (size_t i = 0; i < nx; i++) {x_centers_.push_back((x_faces_[i] + x_faces_[i+1]) / 2.0);}
    for (size_t j = 0; j < ny; j++) {y_centers_.push_back((y_faces_[j] + y_faces_[j+1]) / 2.0);}
    for (size_t k = 0; k < nz; k++) {z_centers_.push_back((z_faces_[k] + z_faces_[k+1]) / 2.0);}
}
const std::vector<double> UniformMesh::x_centers() const {return x_centers_;}
const std::vector<double> UniformMesh::y_centers() const {return y_centers_;}
const std::vector<double> UniformMesh::z_centers() const {return z_centers_;}
const std::vector<double>   UniformMesh::x_faces() const {return   x_faces_;}
const std::vector<double>   UniformMesh::y_faces() const {return   y_faces_;}
const std::vector<double>   UniformMesh::z_faces() const {return   z_faces_;}
size_t UniformMesh::get_nx() const {return x_centers_.size();}
size_t UniformMesh::get_ny() const {return y_centers_.size();}
size_t UniformMesh::get_nz() const {return z_centers_.size();}
double UniformMesh::get_hx() const {return (xlim1_ - xlim0_) / get_nx();}
double UniformMesh::get_hy() const {return (ylim1_ - ylim0_) / get_ny();}
double UniformMesh::get_hz() const {return (zlim1_ - zlim0_) / get_nz();}
size_t UniformMesh::get_ndim() const {
    if      (get_nz() > 1) {return 3;}
    else if (get_ny() > 1) {return 2;}
    else                   {return 1;}
}