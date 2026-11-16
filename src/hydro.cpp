#include <algorithm>
#include "hydro.hpp"

double prim_to_cs(const Primitives& w, double gamma) {
    double cs_sq = gamma * w.P / w.rho;
    return std::sqrt(std::max(cs_sq, 0.0));
}

Primitives cons_to_prim(const Conservatives& u, double gamma) {
    double v1 = u.m1 / u.rho;
    double v2 = u.m2 / u.rho;
    double v3 = u.m3 / u.rho;
    double KE = 0.5 * ((u.m1 * u.m1) + (u.m2 * u.m2) + (u.m3 * u.m3));
    double P = (gamma - 1.0) * (u.E - KE);
    Primitives w{u.rho, v1, v2, v3, P};
    return w;
}

Conservatives prim_to_cons(const Primitives& w, double gamma) {
    double m1 = w.rho * w.v1;
    double m2 = w.rho * w.v2;
    double m3 = w.rho * w.v3;
    double KE = 0.5 * w.rho * ((w.v1 * w.v1) + (w.v2 * w.v2) + (w.v3 * w.v3));
    double E = KE + (w.P / (gamma - 1.0));
    Conservatives u{w.rho, m1, m2, m3, E};
    return u;
}

Conservatives prim_to_flux(const Primitives& w, double gamma) {
    double frho = w.rho * w.v1;
    double fm1 = (w.rho * w.v1 * w.v1) + w.P;
    double fm2 = w.rho * w.v1 * w.v2;
    double fm3 = w.rho * w.v1 * w.v3;
    double KE = 0.5 * w.rho * ((w.v1 * w.v1) + (w.v2 * w.v2) + (w.v3 * w.v3));
    double fE = (KE + (w.P / (gamma - 1.0)) + w.P) * w.v1;
    Conservatives u{frho, fm1, fm2, fm3, fE};
    return u;
}
