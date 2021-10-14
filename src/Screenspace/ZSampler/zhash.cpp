
#include <cstdint>
#include <cassert>
#include <algorithm>

#include "ztools.h"
#include "zcommon.h"

static
uint32_t hash(uint32_t seqNo, uint32_t dim) {
    const int BITS = 24;
    const uint32_t MASK = (1 << BITS) - 1;
    const uint32_t Z = 0x9E377A;                                                // Z / (1 << BITS) approximates 1 - golden ratio
    seqNo ^= dim * 0x555555;                                                    // This should embed the bits of dim with all the bits of seqNo
    uint32_t x = (seqNo * Z) & MASK;                                            // Fractional part
    return (x * 24) >> BITS;                                                    // Map to the desired range
}

static
TileInfo getChildInfo(
    uint32_t X, uint32_t Y,                                                     // Absolute coordinates of child tile.
    TileInfo rootInfo,                                                          // ID and seqNo of root tile
    int digits,                                                                 // Subdivision depth of parent tile; only these digits are considered!
    int dim                                                        // Array of bytes indicating the ranks of children for each ID
) {
    uint32_t tileID = rootInfo.id;
    uint32_t seqNo = rootInfo.seqNo;
    int shift = (digits - 1) * shift1D;
    for (int i = 0; i < digits; i++) {
        seqNo <<= shift2D;
        unsigned y = (Y >> shift) & mask1D;
        unsigned x = (X >> shift) & mask1D;
        unsigned childNo = (y << shift1D) | x;
        uint32_t rank = childNo2rank[ hash(tileID, dim) ][childNo];
        seqNo |= rank;
        tileID = (tileID << shift2D) | childNo;
        shift -= shift1D;
    }
    return {tileID, seqNo};
}

static
uint32_t seqNo2position(
    uint32_t seqNo,
    uint32_t tileID,                                                            // ID of root node
    int bits,                                                                   // depth in bits rather than base-4 digits, to handle odd powers of two
    int dim                                                                     // Sampled dimension
) {
    unsigned position = 0;
    if (bits & 1) {                                                             // Special treatment for an odd power of two, hence not a power of 4
        bits--;
        unsigned rank = (seqNo >> bits) & mask1D;                               // Theoretically, masking should not be needed
        unsigned childNo = rank ^ (rank2childNo[ hash(tileID, dim) ][0] & 1);       // We choose only child 0 or 1, and we use the LSB of rank1 as a permutation
        position |= childNo << bits;
    }
    while (bits > 0) {
        bits -= 2;
        unsigned rank = (seqNo >> bits) & mask2D;
        unsigned childNo = rank2childNo[ hash(tileID, dim) ][rank];
        position |= childNo << bits;
        tileID = (tileID << shift2D) | childNo;
    }
    return position;
}

static
Float fixedPt2Float(uint32_t v) {
    return std::min(v * ONE_IN_2E32, OneMinusEpsilon);
}

void Z2DHash(
    int px, int py,                                                                 // Pixel coordinates
    int nSamplesPerPixelSample,
    int nPixelSamples,
    int depth,
    Float *samples,
    int dim
) {
    TileInfo pixelTile = getChildInfo(                                          // Retrieve ID and seqNo of pixel tile
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
        uint32_t position = seqNo2position(sampleNo, pixelTile.id, bits, dim);
        samples[2*position] = fixedPt2Float(fixedPt1);
        samples[2*position+1] = fixedPt2Float(fixedPt2);
        fixedPt1 ^= ZMatrix1stD[CountTrailingZeros(i + 1)];
        fixedPt2 ^= ZMatrix2ndD[CountTrailingZeros(i + 1)];
    }
}

