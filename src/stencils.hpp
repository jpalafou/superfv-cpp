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
    }
    namespace center {
        inline const Stencil p0_1{1};
        inline const Stencil p2_3{-Real{1} / Real{24}, Real{13} / Real{12}, -Real{1} / Real{24}};
    }
    namespace right {
        inline const Stencil p0(left::p0.weights);
        inline const Stencil p1{std::vector<Real>(left::p1.weights.rbegin(), left::p1.weights.rend())};
        inline const Stencil p2{std::vector<Real>(left::p2.weights.rbegin(), left::p2.weights.rend())};
        inline const Stencil p3{std::vector<Real>(left::p3.weights.rbegin(), left::p3.weights.rend())};
    }
    const Stencil& get(Position pos, size_t p);
}