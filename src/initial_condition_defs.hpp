#pragma once
#include "array.hpp"
#include "mesh.hpp"

enum class Variable {rho=0, vx=1, vy=2, vz=3, P=4};

enum Setup {SINUS, SQUARE};

double sinus(double x, double y, double z, Variable v, size_t ndim, const std::vector<double>& params);
double square(double x, double y, double z, Variable v, size_t ndim, const std::vector<double>& params);
double ic_func(double x, double y, double z, Variable v, size_t ndim, const std::vector<double>& params, Setup setup);