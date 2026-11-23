#pragma once
#include "array.hpp"
#include "common.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"

Array ic(const UniformMesh& mesh, size_t p, Setup setup, Real gamma, const std::vector<Real>& params = {});