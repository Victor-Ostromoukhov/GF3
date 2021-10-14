
#include <cstdint>
#include <cassert>
#include <vector>
#include <algorithm>

#include "../../Samplers/Random.h"
#include "ztools.h"


template <typename T>
void Shuffle( T *samples, int count, int nDimensions, RNG &rng ) 
{
    for (int i = 0; i < count; ++i) 
    {
        int other = i + rng.sample_range(count - i);
        for (int j = 0; j < nDimensions; ++j)
            std::swap(samples[nDimensions * i + j], samples[nDimensions * other + j]);
    }
}

struct Point2u
{
    uint32_t x, y;
};

void GrayCodeSample( const uint32_t *C0, const uint32_t *C1, uint32_t n, const Point2u &scramble, Point2f *p ) 
{
    uint32_t v[2] = { scramble.x, scramble.y };
    
    for (uint32_t i = 0; i < n; ++i) {
        p[i].x = std::min(v[0] * Float(2.3283064365386963e-10), OneMinusEpsilon);
        p[i].y = std::min(v[1] * Float(2.3283064365386963e-10), OneMinusEpsilon);
        v[0] ^= C0[CountTrailingZeros(i + 1)];
        v[1] ^= C1[CountTrailingZeros(i + 1)];
    }
}

void Sobol2D( int nSamplesPerPixelSample, int nPixelSamples, Point2f *samples, RNG &rng ) 
{
    Point2u scramble;
    scramble.x= rng.sample();
    scramble.y= rng.sample();

    // Define 2D Sobol$'$ generator matrices _CSobol[2]_
    const uint32_t CSobol[2][32] = 
    {
        {
            0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000,
            0x2000000, 0x1000000, 0x800000, 0x400000, 0x200000, 0x100000, 0x80000,
            0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000, 0x1000, 0x800,
            0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1
        },
        {
            0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000,
            0xaa000000, 0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000,
            0x88880000, 0xcccc0000, 0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000,
            0xa000a000, 0xf000f000, 0x88008800, 0xcc00cc00, 0xaa00aa00, 0xff00ff00,
            0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0, 0x88888888, 0xcccccccc,
            0xaaaaaaaa, 0xffffffff
        }
    };
    
    GrayCodeSample(CSobol[0], CSobol[1], nSamplesPerPixelSample * nPixelSamples, scramble, samples);
    
    for (int i = 0; i < nPixelSamples; ++i)
        Shuffle(samples + i * nSamplesPerPixelSample, nSamplesPerPixelSample, 1, rng);
    Shuffle(samples, nPixelSamples, nSamplesPerPixelSample, rng);
}

void ZTSequence2D(
    int nSamplesPerPixelSample,
    int nPixelSamples,
    std::vector< std::vector<Point2f> >& samples2D,
    int nDimensions,
    unsigned int seed
) {
    assert(int(samples2D.size()) == nDimensions);
    RNG rng(seed);
    for (size_t i = 0; i < samples2D.size(); ++i)
        Sobol2D(1, nPixelSamples, &samples2D[i][0], rng);
}

