#include <algorithm>
#include "hydro.hpp"

Conservatives llf(Primitives& wl, Primitives& wr, double gamma) {
    double       csl =   prim_to_cs(wl, gamma);
    double       csr =   prim_to_cs(wr, gamma);
    Conservatives ul = prim_to_cons(wl, gamma);
    Conservatives ur = prim_to_cons(wr, gamma);
    Conservatives fl = prim_to_flux(wl, gamma);
    Conservatives fr = prim_to_flux(wr, gamma);
    double sl = csl + std::abs(wl.v1);
    double sr = csr + std::abs(wr.v1);
    double smax = std::max(sl, sr);

    double frho = (fl.rho + fr.rho - (smax * (ur.rho - ul.rho)));
    double  fm1 = (fl.m1  + fr.m1  - (smax * (ur.m1  - ul.m1)));
    double  fm2 = (fl.m2  + fr.m2  - (smax * (ur.m2  - ul.m2)));
    double  fm3 = (fl.m3  + fr.m3  - (smax * (ur.m3  - ul.m3)));
    double   fE = (fl.E   + fr.E   - (smax * (ur.E   - ul.E)));
    Conservatives fluxes{frho, fm1, fm2, fm3, fE};
    return fluxes;
}