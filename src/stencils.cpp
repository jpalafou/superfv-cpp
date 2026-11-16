#include <cstddef>
#include <vector>
#include "stencils.hpp"

const double& Stencil::operator[](std::size_t i) const {return data[i];}
const std::size_t Stencil::size() const {return data.size();}
const std::size_t Stencil::left_reach() const {return data.size() / 2;}
const std::size_t Stencil::right_reach() const {return (data.size() - 1) / 2;}

const Stencil conservative_interpolation(int p, Position pos) {   
    std::vector<double> data; 
    switch (p) {
        case 0:
            data.push_back(1.0);
            break;
        case 1:
            if (pos == Position::left) {
                data.push_back(1.0 / 4.0);
                data.push_back(1.0);
                data.push_back(-1.0 / 4.0);
            } else if (pos == Position::center) {
                data.push_back(1.0);
            } else {
                data.push_back(-1.0 / 4.0);
                data.push_back(1.0);
                data.push_back(1.0 / 4.0);
            }
            break;
        default:
            throw std::invalid_argument("p.");
            
    }
    const Stencil stencil{data, (data.size() - 1) / 2};
    return stencil;
}