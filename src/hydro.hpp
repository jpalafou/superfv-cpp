#pragma once
#include "common.hpp"

struct HydroParams{const Real gamma, cfl;};

struct    Primitives {Real rho, v1, v2, v3, P;};
struct Conservatives {Real rho, m1, m2, m3, E;};

Real            prim_to_cs(const Primitives& w,    Real gamma);
Primitives    cons_to_prim(const Conservatives& u, Real gamma);
Conservatives prim_to_cons(const Primitives& w,    Real gamma);
Conservatives prim_to_flux(const Primitives& w,    Real gamma);