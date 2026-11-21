#pragma once
#include <vector>
#include "common.hpp"

struct Stencil {
    const std::vector<Real> data;
    const std::size_t center_idx;

    const Real& operator[](std::size_t i) const;
    const std::size_t size() const;
    const std::size_t left_reach() const;
    const std::size_t right_reach() const;
};

enum Position {left, center, right};

const Stencil conservative_interpolation(int p, Position pos);