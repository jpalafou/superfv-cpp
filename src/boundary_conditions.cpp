#include <stdexcept>
#include "array.hpp"
#include "boundary_conditions.hpp"

Array get_padded_u(const Array& u, size_t nghost){
    switch (u.get_ndim()) {
        case 2:
            return Array(u.get_nvars(), u.get_nx() + 2 * nghost);
        case 3:
            return Array(u.get_nvars(), u.get_nx() + 2 * nghost, u.get_ny() + 2 * nghost);
        case 4:
            return Array(u.get_nvars(), u.get_nx() + 2 * nghost, u.get_ny() + 2 * nghost, u.get_nz() + 2 * nghost);
        default:
            throw std::invalid_argument("u must have ndim 2, 3, or 4.");
    }
}

void copy_inner_1d(const Array& u, Array& ubc, size_t nghost) {
    for (size_t v = 0; v < u.get_nvars(); v++) {
        for (size_t i = 0; i < u.get_nx(); i++) {
            ubc(v, i + nghost) = u(v, i);
        }
    }
    return;
}

void copy_inner_2d(const Array& u, Array& ubc, size_t nghost) {
    for (size_t v = 0; v < u.get_nvars(); v++) {
        for (size_t i = 0; i < u.get_nx(); i++) {
            for (size_t j = 0; j < u.get_ny(); j++) {
                ubc(v, i + nghost, j + nghost) = u(v, i, j);
            }
        }
    }
    return;
}

void copy_inner_3d(const Array& u, Array& ubc, size_t nghost) {
    for (size_t v = 0; v < u.get_nvars(); v++) {
        for (size_t i = 0; i < u.get_nx(); i++) {
            for (size_t j = 0; j < u.get_ny(); j++) {
                for (size_t k = 0; k < u.get_nz(); k++) {
                    ubc(v, i + nghost, j + nghost) = u(v, i, j, k);
                }
            }
        }
    }
    return;
}

void periodic_bc_applier(const Array& u, Array& ubc, size_t nghost, Direction direction, size_t v, size_t i, size_t j, size_t k) {
    switch (direction) {
        case Direction::x:
            ubc(v, i, j, k) = u(v, u.get_nx() - nghost + i, j, k);
            ubc(v, ubc.get_nx() - nghost + i, j, k) = u(v, i, j, k);
            return;
        case Direction::y:
            ubc(v, i, j, k) = u(v, i, u.get_ny() - nghost + j, k);
            ubc(v, i, ubc.get_ny() - nghost + j, k) = u(v, i, j, k);
            return;
        default:
            ubc(v, i, j, k) = u(v, i, j, u.get_nz() - nghost + k);
            ubc(v, i, j, ubc.get_nz() - nghost + k) = u(v, i, j, k);
    }
}

void bc_applier(const Array& u, Array& ubc, size_t nghost, Direction direction, Boundary boundary, size_t v, size_t i, size_t j, size_t k, BC bc) {
    switch (bc) {
        case BC::none:
            return;
        case BC::periodic:
            periodic_bc_applier(u, ubc, nghost, direction, v, i, j, k);
            return;
        default:
            throw std::invalid_argument("invalid bc type.");
    }
}

// x  x  .  .  .  .  .  .  x  x

void apply_bc_1d(const Array& u, Array& ubc, size_t nghost, BC lower, BC upper) {
    for (size_t v = 0; v < u.get_nvars(); v++) {
        for (size_t i = 0; i < nghost; i++) {
            bc_applier(u, ubc, nghost, Direction::x, Boundary::lower, v, i, 0, 0, lower);
        }
        for (size_t i = ubc.get_nx() - nghost; i < ubc.get_nx(); i++) {
            bc_applier(u, ubc, nghost, Direction::x, Boundary::upper, v, i, 0, 0, upper);
        }
    }
}

// y  y  y  y  y  y  y  y  y  y
// y  y  y  y  y  y  y  y  y  y
// x  x  .  .  .  .  .  .  x  x
// x  x  .  .  .  .  .  .  x  x
// x  x  .  .  .  .  .  .  x  x
// x  x  .  .  .  .  .  .  x  x
// x  x  .  .  .  .  .  .  x  x
// x  x  .  .  .  .  .  .  x  x
// y  y  y  y  y  y  y  y  y  y
// y  y  y  y  y  y  y  y  y  y

void apply_bc_2d(const Array& u, Array& ubc, size_t nghost, BC x_lower, BC x_upper, BC y_lower, BC y_upper) {
    for (size_t v = 0; v < u.get_nvars(); v++) {
        for (size_t j = nghost; j < ubc.get_ny() - nghost; j++) {
            for (size_t i = 0; i < nghost; i++) {
                bc_applier(u, ubc, nghost, Direction::x, Boundary::lower, v, i, j, 0, x_lower);
            }
            for (size_t i = ubc.get_nx() - nghost; i < ubc.get_nx(); i++) {
                bc_applier(u, ubc, nghost, Direction::x, Boundary::upper, v, i, j, 0, x_upper);
            }
        }
        for (size_t i = 0; i < ubc.get_nx(); i++) {
            for (size_t j = 0; j < nghost; j++) {
                bc_applier(u, ubc, nghost, Direction::y, Boundary::lower, v, i, j, 0, y_lower);
            }
            for (size_t j = ubc.get_ny() - nghost; j < ubc.get_ny(); j++) {
                bc_applier(u, ubc, nghost, Direction::y, Boundary::upper, v, i, j, 0, y_upper);
            }
        }
    }
}

void apply_bc(
    const Array& u, 
    Array& ubc, 
    size_t nghost, 
    BC x_lower,
    BC x_upper,
    BC y_lower,
    BC y_upper,
    BC z_lower,
    BC z_upper
) {
    if (!has_same_ndim(u, ubc)) {throw std::invalid_argument("u and ubc must have matching ndim");}
    switch (u.get_ndim()) {
        case 2:
            copy_inner_1d(u, ubc, nghost);
            apply_bc_1d(u, ubc, nghost, x_lower, x_upper);
            return;
        case 3:
            copy_inner_2d(u, ubc, nghost);
            apply_bc_2d(u, ubc, nghost, x_lower, x_upper, y_lower, y_upper);
            return;
        case 4:
            copy_inner_3d(u, ubc, nghost);
            throw std::invalid_argument("apply_bc_3d is not implemented.");
        default:
            throw std::invalid_argument("u must have ndim 2, 3, or 4.");
    }
}