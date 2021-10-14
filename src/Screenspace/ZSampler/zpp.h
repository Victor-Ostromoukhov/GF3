
#ifndef _ZPP_SAMPLER_H
#define _ZPP_SAMPLER_H

#include "ztools.h"

struct ZPPTable
{
    std::vector<unsigned int> production;
    std::vector<unsigned char> allRanks;
    int tileCount;
    int maxdim;
};

void initZPPTable( ZPPTable& table, const int maxdim, unsigned int *production, const int size );
void initZPPTable( ZPPTable& table, const int maxdim= 4 );
void randomZPPTable( ZPPTable& table, const int maxdim= 4 );

void ZPP2D(
    const ZPPTable& table,
    int px, int py,                                                                  // Pixel coordinates
    int nSamplesPerPixelSample,
    int nPixelSamples,
    int depth,
    double *samples,
    int dim
);

#endif


