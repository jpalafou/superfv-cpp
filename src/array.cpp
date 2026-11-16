#include <cstddef>
#include "array.hpp"

Array::Array(size_t n1, size_t n2)                          : nvars(n1), nx(n2), ny(1),  nz(1),  data(n1 * n2) {}
Array::Array(size_t n1, size_t n2, size_t n3)               : nvars(n1), nx(n2), ny(n3), nz(1),  data(n1 * n2 * n3) {}
Array::Array(size_t n1, size_t n2, size_t n3, size_t n4)    : nvars(n1), nx(n2), ny(n3), nz(n4), data(n1 * n2 * n3 * n4) {}

double& Array::operator()(std::size_t v) {return data[v];}
double& Array::operator()(std::size_t v, std::size_t i) {
    return data[(v * nx) + i];
}
double& Array::operator()(std::size_t v, std::size_t i, std::size_t j) {
    return data[(v * nx * ny) + (i * ny) + j];
}
double& Array::operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) {
    return data[(v * nx * ny * nz) + (i * ny * nz) + (j * nz) + k];
}

const double& Array::operator()(std::size_t v) const {return data[v];}
const double& Array::operator()(std::size_t v, std::size_t i) const {
    return data[(v * nx) + i];
}
const double& Array::operator()(std::size_t v, std::size_t i, std::size_t j) const {
    return data[(v * nx * ny) + (i * ny) + j];
}
const double& Array::operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) const {
    return data[(v * nx * ny * nz) + (i * ny * nz) + (j * nz) + k];
}

const std::size_t Array::get_nvars() const {return nvars;}
const std::size_t Array::get_nx() const {return nx;}
const std::size_t Array::get_ny() const {return ny;}
const std::size_t Array::get_nz() const {return nz;}
const std::size_t Array::get_ndim() const {
    if (nz > 1) {return 4;}
    if (ny > 1) {return 3;}
    return 2;
};

bool has_same_shape(const Array& arr1, const Array& arr2) {
    if (arr1.get_nvars() != arr2.get_nvars()) {return false;}
    if (arr1.get_nx() != arr2.get_nx()) {return false;}
    if (arr1.get_ny() != arr2.get_ny()) {return false;}
    if (arr1.get_nz() != arr2.get_nz()) {return false;}
    return true;
}