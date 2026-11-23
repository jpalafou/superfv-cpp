#pragma once
#include <cstddef>
#include <vector>
#include "array.hpp"
#include "common.hpp"
#include "initial_condition_defs.hpp"
#include "mesh.hpp"

Array ic(const UniformMesh& mesh, std::size_t p, Setup setup, Real gamma, const std::vector<Real>& params = {});