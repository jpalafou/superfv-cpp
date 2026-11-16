#pragma once

struct    Primitives {double rho, v1, v2, v3, P;};
struct Conservatives {double rho, m1, m2, m3, E;};

double          prim_to_cs(const Primitives& w,    double gamma);
Primitives    cons_to_prim(const Conservatives& u, double gamma);
Conservatives prim_to_cons(const Primitives& w,    double gamma);
Conservatives prim_to_flux(const Primitives& w,    double gamma);