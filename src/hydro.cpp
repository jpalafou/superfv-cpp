#include <algorithm>
#include <cmath>
#include "common.hpp"
#include "hydro.hpp"

Real prim_to_cs(const Primitives& w, Real gamma) {
    Real cs_sq = gamma * w.P / w.rho;
    return std::sqrt(std::max(cs_sq, Real{0}));
}

Primitives cons_to_prim(const Conservatives& u, Real gamma) {
    Real v1 = u.m1 / u.rho;
    Real v2 = u.m2 / u.rho;
    Real v3 = u.m3 / u.rho;
    Real KE = 0.5 * ((u.m1 * u.m1) + (u.m2 * u.m2) + (u.m3 * u.m3)) / u.rho;
    Real P = (gamma - Real{1}) * (u.E - KE);
    Primitives w{u.rho, v1, v2, v3, P};
    return w;
}

Conservatives prim_to_cons(const Primitives& w, Real gamma) {
    Real m1 = w.rho * w.v1;
    Real m2 = w.rho * w.v2;
    Real m3 = w.rho * w.v3;
    Real KE = 0.5 * w.rho * ((w.v1 * w.v1) + (w.v2 * w.v2) + (w.v3 * w.v3));
    Real E = KE + (w.P / (gamma - Real{1}));
    Conservatives u{w.rho, m1, m2, m3, E};
    return u;
}

Conservatives prim_to_flux(const Primitives& w, Real gamma) {
    Real frho = w.rho * w.v1;
    Real fm1 = (w.rho * w.v1 * w.v1) + w.P;
    Real fm2 = w.rho * w.v1 * w.v2;
    Real fm3 = w.rho * w.v1 * w.v3;
    Real KE = 0.5 * w.rho * ((w.v1 * w.v1) + (w.v2 * w.v2) + (w.v3 * w.v3));
    Real fE = (KE + (w.P / (gamma - Real{1})) + w.P) * w.v1;
    Conservatives u{frho, fm1, fm2, fm3, fE};
    return u;
}
