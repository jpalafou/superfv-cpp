#pragma once
#include "array.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"

Array ic(const UniformMesh& mesh, Setup setup, const std::vector<double>& params);