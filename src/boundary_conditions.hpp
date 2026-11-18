#pragma once
#include "array.hpp"

enum Direction {x, y, z};
enum Boundary {lower, upper};
enum BC {none, periodic};

Array get_padded_u(const Array& u, size_t nghost);

void apply_bc(
    const Array& u, 
    Array& ubc, 
    size_t nghost, 
    BC x_lower,
    BC x_upper,
    BC y_lower = BC::none,
    BC y_upper = BC::none,
    BC z_lower = BC::none,
    BC z_upper = BC::none
);