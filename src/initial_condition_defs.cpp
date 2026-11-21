#include <cmath>
#include <numbers>
#include <stdexcept>
#include "common.hpp"
#include "initial_condition_defs.hpp"

void validate_params(const std::vector<Real>& params, size_t n) {
    if (!params.empty() && params.size() != n) {throw std::invalid_argument("params must be empty or have correct length.");}
}

Real sinus(Real x, Real y, Real z, Variable v, size_t ndim, const std::vector<Real>& params) {
    validate_params(params, 6);
    Real rho_max, rho_min, uniform_vx, uniform_vy, uniform_vz, uniform_P;
    if (params.empty()) {
        rho_max    = 2.0;
        rho_min    = 1.0;
        uniform_vx = 1.0;
        uniform_vy = 0.0;
        uniform_vz = 0.0;
        uniform_P  = 1.0;
    } else {
        rho_max    = params[0];
        rho_min    = params[1];
        uniform_vx = params[2];
        uniform_vy = params[3];
        uniform_vz = params[4];
        uniform_P  = params[5];
    }
    switch (v) {
        case Variable::rho: {
            Real r;
            if      (ndim == 1) {r = x;}
            else if (ndim == 2) {r = x + y;}
            else                {r = x + y + z;}
            return rho_min + (rho_max - rho_min) * (Real{0.5} * std::sin(Real{2} * Real{std::numbers::pi} * r) + Real{0.5});
        }
        case Variable::vx:
            return uniform_vx;
        case Variable::vy:
            return uniform_vy;
        case Variable::vz:
            return uniform_vz;
        case Variable::P:
            return uniform_P;
        default:
            throw std::invalid_argument("invalid Variable for sinus ic.");
    }
}

Real square(Real x, Real y, Real z, Variable v, size_t ndim, const std::vector<Real>& params) {
    validate_params(params, 6);
    Real rho_max, rho_min, uniform_vx, uniform_vy, uniform_vz, uniform_P;
    if (params.empty()) {
        rho_max    = 2.0;
        rho_min    = 1.0;
        uniform_vx = 1.0;
        uniform_vy = 0.0;
        uniform_vz = 0.0;
        uniform_P  = 1.0;
    } else {
        rho_max    = params[0];
        rho_min    = params[1];
        uniform_vx = params[2];
        uniform_vy = params[3];
        uniform_vz = params[4];
        uniform_P  = params[5];
    }
    switch (v) {
        case Variable::rho: {
            bool inside_x =                           x > 0.25 && x < 0.75;
            bool inside_y = (ndim >= 2) ? inside_x && y > 0.25 && y < 0.75 : true;
            bool inside_z = (ndim == 3) ? inside_y && z > 0.25 && z < 0.75 : true;
            if (inside_x && inside_y && inside_z) {
                return rho_max;
            } else {
                return rho_min;
            }
        }
        case Variable::vx:
            return uniform_vx;
        case Variable::vy:
            return uniform_vy;
        case Variable::vz:
            return uniform_vz;
        case Variable::P:
            return uniform_P;
        default:
            throw std::invalid_argument("invalid Variable for sinus ic.");
    }
}

Real ic_func(Real x, Real y, Real z, Variable v, size_t ndim, const std::vector<Real>& params, Setup setup) {
    switch (setup) {
        case Setup::SINUS:
            return  sinus(x, y, z, v, ndim, params);
        case Setup::SQUARE:
            return square(x, y, z, v, ndim, params);
        default:
            throw std::invalid_argument("invalid Setup.");
    }
}