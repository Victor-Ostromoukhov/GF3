
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <random>

#include "z.h"
#include "zcommon.h"


static
TileInfo getChildInfo(
    const ZTable& table,
    uint32_t X, uint32_t Y,                                                     // Absolute coordinates of child tile.
    TileInfo rootInfo,                                                          // ID and seqNo of root tile
    int digits,                                                                 // Subdivision depth of parent tile; only these digits are considered!
    int dim                                                        // Array of bytes indicating the ranks of children for each ID
) {
    const unsigned char *ranks = &table.allRanks[dim * table.tileCount];
    uint32_t tileID = rootInfo.id;
    uint32_t seqNo = rootInfo.seqNo;
    int shift = (digits - 1) * shift1D;
    for (int i = 0; i < digits; i++) {
        seqNo <<= shift2D;
        unsigned y = (Y >> shift) & mask1D;
        unsigned x = (X >> shift) & mask1D;
        unsigned childNo = (y << shift1D) | x;
        uint32_t rank = childNo2rank[ ranks[tileID] ][childNo];
        seqNo |= rank;
        tileID = table.production[tileID * 4 + childNo];
        shift -= shift1D;
    }
    return {tileID, seqNo};
}

static
uint32_t seqNo2position(
    const ZTable& table,
    uint32_t seqNo,
    uint32_t tileID,                                                            // ID of root node
    int bits,                                                                   // depth in bits rather than base-4 digits, to handle odd powers of two
    int dim                                                                     // Sampled dimension
) {
    const unsigned char *ranks = &table.allRanks[dim * table.tileCount];
    unsigned position = 0;
    if (bits & 1) {                                                             // Special treatment for an odd power of two, hence not a power of 4
        bits--;
        unsigned rank = (seqNo >> bits) & mask1D;                               // Theoretically, masking should not be needed
        unsigned childNo = rank ^ (rank2childNo[ ranks[tileID] ][0] & 1);       // We choose only child 0 or 1, and we use the LSB of rank1 as a permutation
        position |= childNo << bits;
    }
    while (bits > 0) {
        bits -= 2;
        unsigned rank = (seqNo >> bits) & mask2D;
        unsigned childNo = rank2childNo[ ranks[tileID] ][rank];
        position |= childNo << bits;
        tileID = table.production[tileID * 4 + childNo];
    }
    return position;
}

static
Float fixedPt2Float(uint32_t v) {
    return std::min(v * ONE_IN_2E32, OneMinusEpsilon);
}

void Z2D(
    const ZTable& table,
    int px, int py,                                                                  // Pixel coordinates
    int nSamplesPerPixelSample,
    int nPixelSamples,
    int depth,
    Float *samples,
    int dim
) {
    TileInfo pixelTile = getChildInfo(                                          // Retrieve ID and seqNo of pixel tile
        table,
        px, py,                                                               // Pixel coordinates
        {0, 1},                                                                 // A root tile representing the whole film
        depth,                                                                  // Subdivide down to pixel level
        dim
    );
    uint32_t total = nPixelSamples * nSamplesPerPixelSample;
    int bits = CountTrailingZeros(total);
    uint32_t base = pixelTile.seqNo << bits;
    uint32_t fixedPt1 = Sobol32(base, ZMatrix1stD);
    uint32_t fixedPt2 = Sobol32(base, ZMatrix2ndD);
    for (uint32_t i = 0; i < total; i++) {
        uint32_t sampleNo = i ^ (i >> 1);                                       // Generate a Gray code
        //~ uint32_t seqNo = base + sampleNo;
        uint32_t position = seqNo2position(table, sampleNo, pixelTile.id, bits, dim);
        samples[2*position] = fixedPt2Float(fixedPt1);
        samples[2*position+1] = fixedPt2Float(fixedPt2);
        fixedPt1 ^= ZMatrix1stD[CountTrailingZeros(i + 1)];
        fixedPt2 ^= ZMatrix2ndD[CountTrailingZeros(i + 1)];
    }
}


void initZTable( ZTable& table, const int maxdim, unsigned int *production, const int size )
{
    table.maxdim= maxdim;
    table.tileCount= size;
    
    std::random_device hwseed;    // utilise le generateur materiel
    // pas deterministe, remplacer par sampler::rng() qui change avec le seed global
    std::default_random_engine seed(hwseed());
    
    if(production)
        table.production.assign(production, production + size*4);
        
    else
    {
        table.production.resize(size * 4);
        for(int i= 0; i < size * 4; i++)
            table.production[i]= seed() % size;
    }
    
    const int DIMENSIONS= 1024;
    table.allRanks.resize(DIMENSIONS * size);
    for(int i= 0; i < DIMENSIONS * size; i++)
        table.allRanks[i]= seed() % 24;
}

extern 
unsigned ART_2x2_PRODUCTION[4096 * 4];

void initZTable( ZTable& table, const int maxdim )
{
    initZTable(table, maxdim, ART_2x2_PRODUCTION, 4096);
}

void randomZTable( ZTable& table, const int maxdim )
{
    initZTable(table, maxdim, nullptr, 4096);
}

