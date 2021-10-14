
#include "sampler_s19.h"


#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp.hpp"
#include "../Screenspace/samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp.hpp"

Sampler19::Sampler19( const std::vector<int>& d, const int width, const int height, const size_t spp, const unsigned int seed ) : Sampler(d.size(), width*height*spp, seed),
    m_dimensions(d), m_optimized(nullptr), m_width(width), m_height(height), m_spp(spp)
{
    assert(d.size() <= 8);
    m_optimized= nullptr;
    switch(spp)
    {
        case   1: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp;   break;
        case   2: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp;   break;
        case   4: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp;   break;
        case   8: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp;   break;
        case  16: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp;  break;
        case  32: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp;  break;
        case  64: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp;  break;
        case 128: m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp; break;
        default:  m_optimized= samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp; break;
    }
}
