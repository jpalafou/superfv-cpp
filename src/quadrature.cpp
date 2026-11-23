#include <stdexcept>
#include "quadrature.hpp"

namespace gauss_legendre {
    Quadrature get(size_t n) {
        switch (n) {
            case 2:
                return q2;
            case 3:
                return q3;
            case 4:
                return q4;
            default:
                throw std::invalid_argument("Gauss-Legendre quadrature size not supported.");
        }
    }
}