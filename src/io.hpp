#pragma once
#include <cstddef>
#include <limits>
#include "array.hpp"
#include "common.hpp"

struct IOParams {const std::size_t fprint;};

struct StepSummary {
    size_t step_number;
    Real t;
    Real rhomax = std::numeric_limits<Real>::lowest();
    Real rhomin = std::numeric_limits<Real>::max();
    Real Etotal = 0.0;
};

StepSummary get_step_summary(const Array& u, size_t step_number, Real t);

void printer(const StepSummary& step_summary, Real Eref);
