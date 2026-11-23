#include "array.hpp"

#include <cstddef>
#include <stdexcept>
#include "common.hpp"

Array::Array()                                              : nvars(0),  nx(0),  ny(0),  nz(0),  data(0) {}
Array::Array(size_t n1)                                     : nvars(n1), nx(1),  ny(1),  nz(1),  data(n1) {}
Array::Array(size_t n1, size_t n2)                          : nvars(n1), nx(n2), ny(1),  nz(1),  data(n1 * n2) {}
Array::Array(size_t n1, size_t n2, size_t n3)               : nvars(n1), nx(n2), ny(n3), nz(1),  data(n1 * n2 * n3) {}
Array::Array(size_t n1, size_t n2, size_t n3, size_t n4)    : nvars(n1), nx(n2), ny(n3), nz(n4), data(n1 * n2 * n3 * n4) {}

Real& Array::operator()(std::size_t v) {return data[v];}
Real& Array::operator()(std::size_t v, std::size_t i) {
    return data[(v * nx) + i];
}
Real& Array::operator()(std::size_t v, std::size_t i, std::size_t j) {
    return data[(v * nx * ny) + (i * ny) + j];
}
Real& Array::operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) {
    return data[(v * nx * ny * nz) + (i * ny * nz) + (j * nz) + k];
}

const Real& Array::operator()(std::size_t v) const {return data[v];}
const Real& Array::operator()(std::size_t v, std::size_t i) const {
    return data[(v * nx) + i];
}
const Real& Array::operator()(std::size_t v, std::size_t i, std::size_t j) const {
    return data[(v * nx * ny) + (i * ny) + j];
}
const Real& Array::operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) const {
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

bool has_same_ndim(const Array& arr1, const Array& arr2) {return arr1.get_ndim() == arr2.get_ndim();}

Array get_empty_copy(const Array& arr) {return Array(arr.get_nvars(), arr.get_nx(), arr.get_ny(), arr.get_nz());}

void copy_contents(const Array& arr1, Array& arr2) {
    if (!has_same_shape(arr1, arr2)) {throw std::invalid_argument("arr1 and arr2 must have same shape.");}
    for (size_t v = 0; v < arr1.get_nvars(); v++) {
        for (size_t i = 0; i < arr1.get_nx(); i++) {
            for (size_t j = 0; j < arr1.get_ny(); j++) {
                for (size_t k = 0; k < arr1.get_nz(); k++) {
                    arr2(v, i, j, k) = arr1(v, i, j, k);
                }
            }
        }
    }
}