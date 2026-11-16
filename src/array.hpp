#pragma once
#include <cstddef>
#include <vector>

class Array {
    public:
        Array(size_t n1);
        Array(size_t n1, size_t n2);
        Array(size_t n1, size_t n2, size_t n3);
        Array(size_t n1, size_t n2, size_t n3, size_t n4);
        double& operator()(std::size_t v);
        double& operator()(std::size_t v, std::size_t i);
        double& operator()(std::size_t v, std::size_t i, std::size_t j);
        double& operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k);
        const double& operator()(std::size_t v) const;
        const double& operator()(std::size_t v, std::size_t i) const;
        const double& operator()(std::size_t v, std::size_t i, std::size_t j) const;
        const double& operator()(std::size_t v, std::size_t i, std::size_t j, std::size_t k) const;
        const std::size_t get_nvars() const;
        const std::size_t get_nx() const;
        const std::size_t get_ny() const;
        const std::size_t get_nz() const;
        const std::size_t get_ndim() const;
    private:
        std::vector<double> data;
        std::size_t nvars, nx, ny, nz;
};

bool has_same_shape(const Array& arr1, const Array& arr2);