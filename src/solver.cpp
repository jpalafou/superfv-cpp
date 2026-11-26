#include "solver.hpp"

#include <cstddef>
#include <vector>
#include "array.hpp"
#include "common.hpp"
#include "hydro_rhs.hpp"
#include "initial_conditions.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "ode.hpp"

namespace {
    void call_printer_when_appropriate(
        const std::vector<StepSummary>& step_log,
        const IOParams& params,
        bool stop
    ) {
        size_t step_number = step_log.back().step_number;
        if (params.fprint > 0 && (step_number % params.fprint == 0 || stop)) {
            printer(step_log.back(), step_log.front());
        }
    }

    void call_snapshot_writer_when_appropriate(
        const Array& u,
        std::vector<Real>& snapshot_schedule,
        size_t& snapshot_idx,
        Real t,
        const SolverParams& params,
        bool stop
    ) {
        if (!snapshot_schedule.empty() && (snapshot_schedule.front() <= t || stop)) {
            while (!snapshot_schedule.empty() && snapshot_schedule.front() <= t) {
                snapshot_schedule.erase(snapshot_schedule.begin());
            }
            write_snapshot(params.io.output_dir, u, params.mesh, t, snapshot_idx);
            snapshot_idx++;
        }
    }
}

Array superfv(const SolverParams& params) {
    UniformMesh mesh = UniformMesh(params.mesh);
    Array uold = ic(mesh, params.p, params.setup, params.hydro.gamma);
    Array f(uold.get_nvars(), uold.get_nx() + 1, uold.get_ny(), uold.get_nz());
    Array g, h;
    if (params.mesh.ndim() >= 2) {
        g = Array(uold.get_nvars(), uold.get_nx(), uold.get_ny() + 1, uold.get_nz());
    } 
    if (params.mesh.ndim() == 3) {
        h = Array(uold.get_nvars(), uold.get_nx(), uold.get_ny(), uold.get_nz() + 1);
    }
    Array unew = get_empty_copy(uold);
    Real t = 0;
    std::size_t step_number = 0;
    std::vector<StepSummary> step_log = {get_step_summary(uold, 0, 0.0)};
    std::size_t snapshot_idx = 0;
    std::vector<Real> snapshot_schedule = get_snapshot_schedule(params.tstop, params.io.nout);

    call_printer_when_appropriate(step_log, params.io, false);
    call_snapshot_writer_when_appropriate(uold, snapshot_schedule, snapshot_idx, 0.0, params, false);

    while (true) {
        bool stop = false;
        
        Real dt = compute_dt(uold, params.hydro.cfl, params.hydro.gamma, params.mesh);
        if (params.clip_last_dt && t + dt > params.tstop) {dt = params.tstop - t; stop=true;}
        
        take_step(uold, unew, f, g, h, dt, params);
        t += dt;
        step_number++;
        step_log.push_back(get_step_summary(unew, step_number, t));
        
        if (t >= params.tstop || step_number >= params.max_steps) {stop=true;}
        call_printer_when_appropriate(step_log, params.io, stop);
        call_snapshot_writer_when_appropriate(unew, snapshot_schedule, snapshot_idx, t, params, stop);
        if (stop) {break;}
        
        copy_contents(unew, uold);
    }

    return unew;
}
