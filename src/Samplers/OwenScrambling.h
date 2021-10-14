//
// Created by lpaulin on 17/09/20.
//

#ifndef SOBOL_OWENSCRAMBLING_H
#define SOBOL_OWENSCRAMBLING_H

#include <bitset>
#include <iostream>
#include "Random.h"

//Default tree depth
#define OWENPLUS_TREE_DEPTH 32u
//OwenPlus Scrambled points with depth of 32 range between 0 and 4294967296
#define UINT32SOBOLNORM 4294967296U

//seed has to stay the same for a given dimention
template <typename integer>
integer OwenScrambling(integer sobolPoint, uint32_t seed, const uint32_t owen_tree_depth) {

    RNG scramble_rng;

    constexpr uint32_t bits= 8*sizeof(integer);
    constexpr uint32_t shift= 8*sizeof(integer) -1;

    // flip root, node_index == 0, implicit rng.index(0)
    scramble_rng.seed(seed);
    integer flip= scramble_rng.sample_range(2) << shift;
    integer code= sobolPoint ^ flip;        // flip MSB

    for(uint32_t idigit= 1; idigit < owen_tree_depth; idigit++)
    {
        integer level_base= (1u << idigit) - 1;
        integer level_offset = sobolPoint >> (bits - idigit);	// level > 0 == 2^d nodes
//if ( idigit == 1 ) std::cout << level_offset << " lobset" << std::endl;
        integer node_index = level_base + level_offset;

        integer node_seed= seed ^ node_index;
        scramble_rng.seed(node_seed);

        integer flip= scramble_rng.sample_range(2) << (shift - idigit);
        code= code ^ flip;
    }

    return code;
}


#endif //SOBOL_OWENSCRAMBLING_H
