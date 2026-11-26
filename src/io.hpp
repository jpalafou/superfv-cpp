#pragma once
#include <cstddef>
#include <filesystem>
#include <limits>
#include <vector>
#include "array.hpp"
#include "common.hpp"
#include "mesh.hpp"

struct IOParams {
    const std::size_t fprint;
    const std::size_t nout;
    const std::filesystem::path output_dir;
};

struct StepSummary {
    size_t step_number;
    Real t;
    Real rhomax = std::numeric_limits<Real>::lowest();
    Real rhomin = std::numeric_limits<Real>::max();
    Real Etotal = 0.0;
};

StepSummary get_step_summary(const Array& u, size_t step_number, Real t);

void printer(const StepSummary& step_summary, Real Eref);

std::vector<Real> get_snapshot_schedule(Real tstop, std::size_t nout);

// Binary layout (native byte order): nx, ny, nz as uint64_t; then
// xlim0, xlim1, ylim0, ylim1, zlim0, zlim1, t, and u as Real.
// u is stored in (variable, x, y, z) order.
void write_snapshot(const std::filesystem::path& output_dir, const Array& u, const UniformMeshParams& mesh, Real t, std::size_t n);
