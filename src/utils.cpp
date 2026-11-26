#include "utils.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>
#include "common.hpp"

std::vector<Real> linspace(Real lower, Real upper, std::size_t n){
    std::vector<Real> out;
    if (n < 2) {throw std::invalid_argument("n must be greater than 1.");}
    Real h = (upper - lower) / (n - 1);
    for(std::size_t i = 0; i < n - 1; ++i) {
        out.push_back(lower + h * i);
    }
    out.push_back(upper);
    return out;
}
