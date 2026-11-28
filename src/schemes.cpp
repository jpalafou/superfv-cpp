#include "schemes.hpp"

#include <stdexcept>

namespace {
    std::vector<BaseScheme> make_fallback_schemes(std::size_t p, FallbackCascade cascade, SlopeLimiter limiter) {
        if (p < 1 && cascade != FallbackCascade::none) {throw std::invalid_argument("need p > 0 for fallback cascade.");}
        switch (cascade) {
            case FallbackCascade::none:
                return {};
            case FallbackCascade::full: {
                std::vector<BaseScheme> schemes;
                for (std::size_t degree = p; degree > 0; --degree) {schemes.emplace_back(degree - 1);}
                return schemes;
            }
            case FallbackCascade::MUSCL:
                return {BaseScheme{limiter}};
            case FallbackCascade::MUSCL0:
                return {BaseScheme{limiter}, BaseScheme{0}};
            default:
                throw std::invalid_argument("unknown fallback cascade");
        }
    }
}

BaseScheme::BaseScheme(  std::size_t degree) : p(degree), MUSCL_limiter(SlopeLimiter::none) {}
BaseScheme::BaseScheme(SlopeLimiter limiter) : p(     1), MUSCL_limiter(           limiter) {}
bool BaseScheme::is_MUSCL() const {return p == 1 && MUSCL_limiter != SlopeLimiter::none;}

MOOD_Scheme::MOOD_Scheme(
    std::size_t base_degree, 
    std::size_t max_revs, 
    FallbackCascade cascade, 
    SlopeLimiter limiter,
    Real tolerance,
    bool use_NAD,
    bool use_PAD,
    bool use_SED
) :
    BaseScheme(base_degree),
    nrev(max_revs),
    fallback_schemes(max_revs > 0 ? make_fallback_schemes(base_degree, cascade, limiter) : std::vector<BaseScheme>{}),
    tol(tolerance),
    numerical_admissibility_detection(use_NAD),
    physical_admissibility_detection( use_PAD),
    smooth_extrema_detection(         use_SED)
{}

bool MOOD_Scheme::is_MUSCL() const {return BaseScheme::is_MUSCL() && nrev == 0;}
