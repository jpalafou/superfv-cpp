#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include "common.hpp"
#include "stencils.hpp"

Stencil::Stencil(          std::vector<Real> weights_vector) : weights(weights_vector) {};
Stencil::Stencil(std::initializer_list<Real>         values) :         weights(values) {};
const Real& Stencil::operator[](size_t i) const {return weights[i];}
size_t        Stencil::size() const {return   weights.size();}
size_t  Stencil::left_reach() const {return       size() / 2;}
size_t Stencil::right_reach() const {return (size() - 1) / 2;}

namespace conservative_interpolation{
    const Stencil& get(Position pos, size_t p) {
        switch (pos) {
            case Position::left: {
                switch (p) {
                    case 0:
                        return left::p0;
                    case 1:
                        return left::p1;
                    case 2:
                        return left::p2;
                    case 3:
                        return left::p3;
                    default:
                        throw std::invalid_argument("unsupported p for left conservative interpolation stencil.");
                }
            }
            case Position::center: {
                switch (p) {
                    case 0:
                    case 1:
                        return center::p0_1;
                    case 2:
                    case 3:
                        return center::p2_3;
                    default:
                        throw std::invalid_argument("unsupported p for center conservative interpolation stencil.");
                }
            }
            case Position::right: {
                switch (p) {
                    case 0:
                        return right::p0;
                    case 1:
                        return right::p1;
                    case 2:
                        return right::p2;
                    case 3:
                        return right::p3;
                    default:
                        throw std::invalid_argument("unsupported p for right conservative interpolation stencil.");
                }
            }
            default:
                throw std::invalid_argument("unsupported position for conservative interpolation.");
        }
    }
}