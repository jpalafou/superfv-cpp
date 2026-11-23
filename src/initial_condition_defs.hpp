#pragma once
#include <cstddef>
#include <vector>
#include "common.hpp"

enum class Variable {rho=0, vx=1, vy=2, vz=3, P=4};

enum Setup {SINUS, SQUARE};

Real sinus(Real x, Real y, Real z, Variable v, std::size_t ndim, const std::vector<Real>& params);
Real square(Real x, Real y, Real z, Variable v, std::size_t ndim, const std::vector<Real>& params);
Real ic_func(Real x, Real y, Real z, Variable v, std::size_t ndim, const std::vector<Real>& params, Setup setup);