#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "common.hpp"

struct Quadrature {const std::vector<Real> x, w;};

namespace gauss_legendre {
    inline const Quadrature q2{
        {
            -std::sqrt(Real{3}) / Real{3},
             std::sqrt(Real{3}) / Real{3}
        },
        {
            Real{1},
            Real{1}
        }
    };
    inline const Quadrature q3{
        {
            -std::sqrt(Real{15}) / Real{5},
             Real{0},
             std::sqrt(Real{15}) / Real{5}
        },
        {
            Real{5} / Real{9},
            Real{8} / Real{9},
            Real{5} / Real{9}
        }
    };
    inline const Quadrature q4{
        {
            -std::sqrt(Real{525} + (Real{70} * std::sqrt(Real{30}))) / Real{35},
            -std::sqrt(Real{525} - (Real{70} * std::sqrt(Real{30}))) / Real{35},
             std::sqrt(Real{525} - (Real{70} * std::sqrt(Real{30}))) / Real{35},
             std::sqrt(Real{525} + (Real{70} * std::sqrt(Real{30}))) / Real{35}
        },
        {
            Real{18 - std::sqrt(Real{30})} / Real{36},
            Real{18 + std::sqrt(Real{30})} / Real{36},
            Real{18 + std::sqrt(Real{30})} / Real{36},
            Real{18 - std::sqrt(Real{30})} / Real{36}
        }
    };
    const Quadrature& get(size_t n);
}