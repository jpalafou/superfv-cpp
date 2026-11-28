#pragma once

#include <cstddef>
#include <vector>
#include "common.hpp"


enum class SlopeLimiter {none, minmod, moncen};
enum class FallbackCascade {none, full, MUSCL, MUSCL0};

struct BaseScheme {
    const std::size_t p;
    const SlopeLimiter MUSCL_limiter;
    
    BaseScheme(std::size_t degree);
    BaseScheme(SlopeLimiter limiter);
    bool is_MUSCL() const;
};

struct MOOD_Scheme : BaseScheme {
    using BaseScheme::BaseScheme;

    const std::size_t nrev = 0;
    const std::vector<BaseScheme> fallback_schemes = {};
    const Real tol = 0.0;
    const bool numerical_admissibility_detection = false;
    const bool physical_admissibility_detection  = false;
    const bool smooth_extrema_detection          = false;
    MOOD_Scheme(
        std::size_t base_degree, 
        std::size_t max_revs, 
        FallbackCascade cascade, 
        SlopeLimiter limiter,
        Real tolerance,
        bool use_NAD,
        bool use_PAD,
        bool use_SED
    );
    bool is_MUSCL() const;
};

