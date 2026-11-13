#include <cstddef>

#include "array.hpp"

Array::Array(size_t n1, size_t n2)                          : nvars(n1), nx(n2), ny(1),  nz(1),  data(n1 * n2) {}
Array::Array(size_t n1, size_t n2, size_t n3)               : nvars(n1), nx(n2), ny(n3), nz(1),  data(n1 * n2 * n3) {}
Array::Array(size_t n1, size_t n2, size_t n3, size_t n4)    : nvars(n1), nx(n2), ny(n3), nz(n4), data(n1 * n2 * n3 * n4) {}

double& Array::operator()(std::size_t v, std::size_t i) {
    return data[(v * nx) + i];
}

double& Array::operator()(std::size_t v, std::size_t i, std::size_t j) {
    return data[(v * nx * ny) + (i * ny) + j];
}

double& Array::operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) {
    return data[(v * nx * ny * nz) + (i * ny * nz) + (j * nz) + k];
}

const double& Array::operator()(std::size_t v, std::size_t i) const {
    return data[(v * nx) + i];
}

const double& Array::operator()(std::size_t v, std::size_t i, std::size_t j) const {
    return data[(v * nx * ny) + (i * ny) + j];
}

const double& Array::operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) const {
    return data[(v * nx * ny * nz) + (i * ny * nz) + (j * nz) + k];
}