
#ifndef _ZSAMPLER_HASH_H
#define _ZSAMPLER_HASH_H

#include "ztools.h"

void Z2DHash(
    int px, int py,                                                                 // Pixel coordinates
    int nSamplesPerPixelSample,
    int nPixelSamples,
    int depth,
    double *samples,
    int dim
);

#endif

