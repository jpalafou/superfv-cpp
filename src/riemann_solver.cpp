#include "riemann_solver.hpp"

#include <algorithm>
#include <cmath>
#include "common.hpp"
#include "hydro.hpp"

Conservatives llf(Primitives& wl, Primitives& wr, Real gamma) {
    Real         csl =   prim_to_cs(wl, gamma);
    Real         csr =   prim_to_cs(wr, gamma);
    Conservatives ul = prim_to_cons(wl, gamma);
    Conservatives ur = prim_to_cons(wr, gamma);
    Conservatives fl = prim_to_flux(wl, gamma);
    Conservatives fr = prim_to_flux(wr, gamma);
    Real sl = csl + std::abs(wl.v1);
    Real sr = csr + std::abs(wr.v1);
    Real smax = std::max(sl, sr);

    Real frho = 0.5 * (fl.rho + fr.rho - (smax * (ur.rho - ul.rho)));
    Real  fm1 = 0.5 * (fl.m1  + fr.m1  - (smax * (ur.m1  - ul.m1)));
    Real  fm2 = 0.5 * (fl.m2  + fr.m2  - (smax * (ur.m2  - ul.m2)));
    Real  fm3 = 0.5 * (fl.m3  + fr.m3  - (smax * (ur.m3  - ul.m3)));
    Real   fE = 0.5 * (fl.E   + fr.E   - (smax * (ur.E   - ul.E)));
    Conservatives fluxes{frho, fm1, fm2, fm3, fE};
    return fluxes;
}