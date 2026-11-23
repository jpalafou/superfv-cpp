#pragma once
#include <cstddef>
#include "array.hpp"

enum Direction {x, y, z};
enum Boundary {lower, upper};
enum BC {none, periodic};

struct BC_Params{
    const std::size_t nghost;
    const BC x_lower;
    const BC x_upper;
    const BC y_lower = BC::none;
    const BC y_upper = BC::none;
    const BC z_lower = BC::none;
    const BC z_upper = BC::none;
};

Array get_padded_u(const Array& u, std::size_t nghost);

void apply_bc(const Array& u, Array& ubc, const BC_Params& params);