
#ifndef SOBOL_CASCADE_GENERATOR_H
#define SOBOL_CASCADE_GENERATOR_H

#include <vector>

#include "SobolGenerator1D.h"
#include "OwenScrambling.h"



template< typename T >
inline void getUniformND( T *values,
    const std::vector<SobolGenerator1D<uint32_t> >& sobols,
    const uint32_t *seeds,
    const int nDims,
    const uint32_t n,
    const uint32_t nbits,
    const uint32_t owen_tree_depth = 32,
    const bool owen_permut_flag = true )
{
    uint32_t IDcode= n;		// dim 0: take n-th point as index -> into 32-bit integer IDcode
    for(int idim = 0; idim < nDims; idim++) 
    {
        IDcode= sobols[idim].getSobolInt(IDcode);	// radix-inversion + sobol
        uint32_t res_IDcode= IDcode;				// used to calculate the resulting value
        IDcode= IDcode >> (32-nbits);				// this will be used as new index for the next dimension
        if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
            res_IDcode= OwenScrambling(res_IDcode, seeds[idim], owen_tree_depth);
        
        values[idim]= T(res_IDcode) / T(UINT32SOBOLNORM);	// final value (double) for this dimension
    }
}


#endif
