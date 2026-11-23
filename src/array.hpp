#pragma once
#include <cstddef>
#include <vector>
#include "common.hpp"

class Array {
    public:
        Array();
        Array(std::size_t n1);
        Array(std::size_t n1, std::size_t n2);
        Array(std::size_t n1, std::size_t n2, std::size_t n3);
        Array(std::size_t n1, std::size_t n2, std::size_t n3, std::size_t n4);
        Real& operator()(std::size_t v);
        Real& operator()(std::size_t v, std::size_t i);
        Real& operator()(std::size_t v, std::size_t i, std::size_t j);
        Real& operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k);
        const Real& operator()(std::size_t v) const;
        const Real& operator()(std::size_t v, std::size_t i) const;
        const Real& operator()(std::size_t v, std::size_t i, std::size_t j) const;
        const Real& operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) const;
        const std::size_t get_nvars() const;
        const std::size_t get_nx() const;
        const std::size_t get_ny() const;
        const std::size_t get_nz() const;
        const std::size_t get_ndim() const;
    private:
        std::vector<Real> data;
        std::size_t nvars, nx, ny, nz;
};

bool has_same_shape(const Array& arr1, const Array& arr2);

bool has_same_ndim(const Array& arr1, const Array& arr2);

Array get_empty_copy(const Array& arr);

void copy_contents(const Array& arr1, Array& arr2);