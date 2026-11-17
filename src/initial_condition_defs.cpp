#include <cmath>
#include <numbers>
#include <stdexcept>
#include "initial_condition_defs.hpp"

void validate_params(const std::vector<double>& params, size_t n) {
    if (!params.empty() && params.size() != n) {throw std::invalid_argument("params must be empty or have correct length.");}
}

double sinus(double x, double y, double z, Variable v, size_t ndim, const std::vector<double>& params) {
    validate_params(params, 6);
    switch (v) {
        case Variable::rho: {
            double rho_max = params[0];
            double rho_min = params[1];
            double r;
            if      (ndim == 1) {r = x;}
            else if (ndim == 2) {r = x + y;}
            else                {r = x + y + z;}
            return rho_min + (rho_max - rho_min) * (0.5 * std::sin(2.0 * std::numbers::pi * r) + 0.5);
        }
        case Variable::vx:
            return params[2];
        case Variable::vy:
            return params[3];
        case Variable::vz:
            return params[4];
        case Variable::P:
            return params[5];
        default:
            throw std::invalid_argument("invalid Variable for sinus ic.");
    }
}

double square(double x, double y, double z, Variable v, size_t ndim, const std::vector<double>& params) {
    validate_params(params, 6);
    switch (v) {
        case Variable::rho: {
            double rho_max = params[0];
            double rho_min = params[1];
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
            return params[2];
        case Variable::vy:
            return params[3];
        case Variable::vz:
            return params[4];
        case Variable::P:
            return params[5];
        default:
            throw std::invalid_argument("invalid Variable for sinus ic.");
    }
}

double ic_func(double x, double y, double z, Variable v, size_t ndim, const std::vector<double>& params, Setup setup) {
    switch (setup) {
        case Setup::SINUS:
            return  sinus(x, y, z, v, ndim, params);
        case Setup::SQUARE:
            return square(x, y, z, v, ndim, params);
        default:
            throw std::invalid_argument("invalid Setup.");
    }
}