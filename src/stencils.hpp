#pragma once
#include <initializer_list>
#include <vector>
#include "common.hpp"

struct Stencil{
    const std::vector<Real> weights;

    Stencil(std::vector<Real> weights_vector);
    Stencil(std::initializer_list<Real> values);
    const Real& operator[](std::size_t i) const;
    size_t        size() const;
    size_t  left_reach() const;
    size_t right_reach() const;
};

enum Position {left, center, right};

namespace conservative_interpolation{
    namespace left {
        inline const Stencil p0{1};
        inline const Stencil p1{Real{1} / Real{4}, 1, -Real{1} / Real{4}};
        inline const Stencil p2{Real{1} / Real{3}, Real{5} / Real{6}, -Real{1} / Real{6}};
        inline const Stencil p3{
            -Real{1} / Real{24},
             Real{5} / Real{12},
             Real{5} /  Real{6},
            -Real{1} /  Real{4},
             Real{1} / Real{24},
        };
        inline const Stencil p4{
             -Real{1} / Real{20},
              Real{9} / Real{20},
             Real{47} / Real{60},
            -Real{13} / Real{60},
              Real{1} / Real{30},
        };
        inline const Stencil p5{
              Real{1} / Real{120},
             -Real{1} /  Real{12},
             Real{59} / Real{120},
             Real{47} /  Real{60},
            -Real{31} / Real{120},
              Real{1} /  Real{15},
             -Real{1} / Real{120},
        };
        inline const Stencil p6{
               Real{1} / Real{105},
             -Real{19} / Real{210},
             Real{107} / Real{210},
             Real{319} / Real{420},
            -Real{101} / Real{420},
               Real{5} /  Real{84},
              -Real{1} / Real{140},
        };
        inline const Stencil p7{
              -Real{1} / Real{560},
              Real{17} / Real{840},
             -Real{97} / Real{840},
             Real{449} / Real{840},
             Real{319} / Real{420},
            -Real{223} / Real{840},
              Real{71} / Real{840},
              -Real{1} /  Real{56},
               Real{1} / Real{560},
        };
    }
    namespace center {
        inline const Stencil p0_1{1};
        inline const Stencil p2_3{-Real{1} / Real{24}, Real{13} / Real{12}, -Real{1} / Real{24}};
        inline const Stencil p4_5{
               Real{3} / Real{640},
             -Real{29} / Real{480},
            Real{1067} / Real{960},
             -Real{29} / Real{480},
               Real{3} / Real{640},
        };
        inline const Stencil p6_7{
               -Real{5} /   Real{7168},
              Real{159} /  Real{17920},
            -Real{7621} / Real{107520},
            Real{30251} /  Real{26880},
            -Real{7621} / Real{107520},
              Real{159} /  Real{17920},
               -Real{5} /   Real{7168},
        };
    }
    namespace right {
        inline const Stencil p0(left::p0.weights);
        inline const Stencil p1{std::vector<Real>(left::p1.weights.rbegin(), left::p1.weights.rend())};
        inline const Stencil p2{std::vector<Real>(left::p2.weights.rbegin(), left::p2.weights.rend())};
        inline const Stencil p3{std::vector<Real>(left::p3.weights.rbegin(), left::p3.weights.rend())};
        inline const Stencil p4{std::vector<Real>(left::p4.weights.rbegin(), left::p4.weights.rend())};
        inline const Stencil p5{std::vector<Real>(left::p5.weights.rbegin(), left::p5.weights.rend())};
        inline const Stencil p6{std::vector<Real>(left::p6.weights.rbegin(), left::p6.weights.rend())};
        inline const Stencil p7{std::vector<Real>(left::p7.weights.rbegin(), left::p7.weights.rend())};
    }
    const Stencil& get(Position pos, size_t p);
}
