
#ifndef _Z_H
#define _Z_H

#include "ztools.h"

struct ZTable
{
    std::vector<unsigned int> production;
    std::vector<unsigned char> allRanks;
    int tileCount;
    int maxdim;
};

void initZTable( ZTable& table, const int maxdim, unsigned int *production, const int size );
void initZTable( ZTable& table, const int maxdim= 4 );
void randomZTable( ZTable& table, const int maxdim= 4 );

void Z2D(
    const ZTable& table,
    int px, int py,                                                                  // Pixel coordinates
    int nSamplesPerPixelSample,
    int nPixelSamples,
    int depth,
    double *samples,
    int dim
);

#endif


