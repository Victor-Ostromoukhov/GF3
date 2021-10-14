#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"


#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <omp.h>

// precalul des rangs des fonctions de tests
//~ #define RANK
#undef RANK

// selectionne les fonctions de tests
#define GAUSS
//~ #undef GAUSS
// gauss ou heaviside


#ifdef RANK
// requires eigen3
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <eigen3/Eigen/Dense>
#endif

#include <CLI11.hpp>
#include "../Integration/Integration.h"

#include "../Samplers/Random.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

#include "../Screenspace/ZSampler/z.h"


double generalizedL2( const std::vector< VecXDynamic > &points, const size_t dimension )
{
    long double a, factor_b, factor_c;
    auto N= points.size();
    auto D= dimension;

    a = pow((4.0/3.0), D);
    factor_b = 2.0/(double)N;

    factor_c = 1.0/(double)(N);
    factor_c *= factor_c;

    long double sumb = 0.0;
#pragma omp parallel for  reduction(+:sumb)
    for(unsigned int i=0; i<N; i++)
    {
        long double prodb = 1.0;
        for(unsigned int j=0; j<D; j++)
        {
            double uij = points[i][j];
            prodb *= ((3.0 - uij*uij)/2.0);
        }
        sumb += prodb;
    }

    long double sumc = 0.0;
#pragma omp parallel for  reduction(+:sumc)
    for(uint i=0; i<N; i++)
    for(unsigned int iprime=0; iprime<N; iprime++)
    {
        long double prodc = 1.0;
        for(unsigned int j=0; j<D; j++)
        {
            double uij = points[i][j];
            double uiprimej = points[iprime][j];
            double m = uij > uiprimej ? uij : uiprimej;//std::max(uij, uiprimej);
            prodc *= (2.0 - m);
        }
        sumc += prodc;
    }

    long double tmp0 = factor_b*sumb;
    long double tmp1 = factor_c*sumc;
    return sqrtl(a -tmp0 + tmp1);
}


// sobol + owenpp / 32 bits 
bool sample_owenpp( std::vector< VecXDynamic >& points, const int spp, const std::vector<int>& dimensions, 
    std::vector< SobolGenerator1D<uint32_t> >& sobols, const int depth )
{
    const int ndim= int(dimensions.size());
    assert(!sobols.empty());
    points.clear();
    VecXDynamic sample(ndim);

    // permutations aleatoires par dimension pour owen
    assert(ndim <= sobols.size());
    std::random_device hwseed;
    for(int d= 0; d < ndim; d++)
        sobols[d].owen_seeds_per_dim= hwseed();
    
    //~ int owen_depth= std::min(int(std::log2(spp)) + depth, 32);
    //~ printf("[owen++] tree depth %d= %d+%d\n", owen_depth, int(std::log2(spp)), depth);
    //~ int owen_depth= std::max(int(std::log2(spp)), depth);
    //~ printf("[owen++] tree depth %d= max(%d, %d)\n", owen_depth, int(std::log2(spp)), depth);
    const int owen_depth= depth;
    
    for(int i= 0; i < spp; i++)
    {
        for(int d= 0; d < ndim; d++)
            sample[d]= getOwenPlus1D(sobols, dimensions[d], i, owen_depth);
            //~ sample[d]= getOwenPlus1D(sobols, d, spp + i, owen_depth); // sobol offset !!
        points.push_back(sample);
    }
    
    return true;
}

// sobol 01 + owenpp / 32 bits 
// n'utilise que les dimensions 01 + permutation differente par paire de dimensions
bool sample_owenpp01( std::vector< VecXDynamic >& points, const int spp, const int ndim, std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    VecXDynamic sample(ndim);

    // permutations aleatoires par dimension pour owen
    assert(ndim <= sobols.size());
    std::random_device hwseed;
    for(int d= 0; d < ndim; d++)
        sobols[d].owen_seeds_per_dim= hwseed();
    
    const int owen_depth= 32;
    // const int owen_depth= std::log2(spp);
    
    for(int i= 0; i < spp; i++)
    {
        for(int pair= 0; pair < ndim / 2; pair++)
        {
            sample[2*pair]= getOwenPlus1D_with_seed(sobols, 0, i, sobols[2*pair].owen_seeds_per_dim, owen_depth);
            sample[2*pair+1]= getOwenPlus1D_with_seed(sobols, 1, i, sobols[2*pair+1].owen_seeds_per_dim, owen_depth);
        }
        
        points.push_back(sample);
    }
    
    return true;
}

// sobol + std owen 
bool sample_owen( std::vector< VecXDynamic >& points, const int spp, const std::vector<int>& dimensions, std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    const int ndim= int(dimensions.size());
    assert(!sobols.empty());
    points.clear();
    VecXDynamic sample(ndim);
    
    // permutations aleatoires par dimension pour owen
    assert(ndim <= sobols.size());
    std::random_device hwseed;
    for(int d= 0; d < ndim; d++)
        sobols[d].owen_seeds_per_dim= hwseed();
    
    int owen_depth= std::log2(spp);
    for(int i= 0; i < spp; i++)
    {
        for(int d= 0; d < ndim; d++)
            sample[d]= getOwenPlus1D(sobols, dimensions[d], i, owen_depth);
        points.push_back(sample);
    }
    
    return true;
}

// sobol 
bool sample_sobol( std::vector< VecXDynamic >& points, const int spp, const std::vector<int>& dimensions, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    const int ndim= int(dimensions.size());
    assert(!sobols.empty());
    points.clear();
    VecXDynamic sample(ndim);
    
    for(int i= 0; i < spp; i++)
    {
        for(int d= 0; d < ndim; d++)
            sample[d]= double(sobols[dimensions[d]].getSobolInt(i)) / double(UINT32SOBOLNORM);
        points.push_back(sample);
    }
    
    return true;
}

// sobol offset
bool sample_sobol_offset( std::vector< VecXDynamic >& points, const int spp, const int ndim, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    VecXDynamic sample(ndim);
    
    for(int i= 0; i < spp; i++)
    {
        for(int d= 0; d < ndim; d++)
            sample[d]= double(sobols[d].getSobolInt(spp+i)) / double(UINT32SOBOLNORM);
        points.push_back(sample);
    }
    
    return true;
}

// cascaded sobol
static void getUniformND(
		VecXDynamic& sample,
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
        const uint32_t nDims,
        const uint32_t n,
		const uint32_t nbits,
        const uint32_t owen_tree_depth = 32,
        const bool owen_permut_flag = true
    ) {
    uint32_t IDcode = n;		// dim 0: take n-th point as index -> into 32-bit integer IDcode
	for (unsigned int idim = 0; idim < nDims; idim++) {
		IDcode = sobols[idim].getSobolInt(IDcode);	// radix-inversion + sobol
		uint32_t res_IDcode = IDcode;				// used to calculate the resulting value
	    IDcode = IDcode >> (32-nbits);				// this will be used as new index for the next dimension
	    if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
        	res_IDcode = OwenScrambling(res_IDcode, sobols[idim].owen_seeds_per_dim, owen_tree_depth);
	    sample[idim] = ((double) res_IDcode / (double) UINT32SOBOLNORM);	// final value (double) for this dimension
	}
}	// getUniformND

bool sample_cascaded_sobol( std::vector< VecXDynamic >& points, const int spp, const int ndim, std::vector< SobolGenerator1D<uint32_t> >& sobols, const int depth )
{
    assert(!sobols.empty());
    points.clear();
    VecXDynamic sample(ndim);
    
    // permutations aleatoires par dimension pour owen
    assert(ndim <= sobols.size());
    std::random_device hwseed;
    for(int d= 0; d < ndim; d++)
        sobols[d].owen_seeds_per_dim= hwseed();
    
    int sobol_bits= std::log2(spp);
    int owen_bits= std::max(sobol_bits + depth, 32);
    
    for(int i= 0; i < spp; i++)
    {
        getUniformND(sample, sobols, ndim, i, sobol_bits, owen_bits, true);
        points.push_back(sample);
    }
    
    return true;
}


// pbrt Z Sampler
bool sample_zsampler( std::vector< VecXDynamic >& points, const int n, const int ndim )
{
    points.clear();
    // verifie que n est un carre 
    int size= std::sqrt(n);
    if(size*size != n)
        return false;   // pas possible
    
    // prealloue tous les samples
    VecXDynamic p(ndim);
    for(int i= 0; i < size*size; i++)
        points.push_back(p);
    
    ZTable ztable;                      // genere une nouvelle table de permutation pour chaque realisation
    initZTable(ztable, ndim/2);         // Zsampler genere des paires de dimensions
    
    int i= 0;
    for(unsigned int y= 0; y < size; y++)
    for(unsigned int x= 0; x < size; x++, i++)
    {
        for(int d= 0; d < ndim / 2; d++)
        {
            double tmp[2];
            Z2D(ztable,
                x, y, // pixel
                1, 1,                   // 1 spp
                Log2Int(size),          // 1 pixel dans une image size x size
                tmp,                    // 1 sample pour la paire de dimensions d
                d);
            
            // re-ordonne les samples
            points[i][2*d]= tmp[0];
            points[i][2*d+1]= tmp[1];
        }
        
        //~ // force la stratification dans le pixel x, y sur une image size x size
        //~ points[i][0]= (x + points[i][0]) / double(size);
        //~ points[i][1]= (y + points[i][1]) / double(size);
    }
    
    assert(points.size() == size*size);
    return true;
}

bool sample_zsampler( std::vector< VecXDynamic >& points, const int size, const int n, const int ndim )
{
    points.clear();
    int nspp= n / (size * size);
    if(nspp != RoundUpPow2(nspp))
        return false;
    if(nspp*size*size != n)
        return false;
    
    // prealloue tous les samples
    VecXDynamic p(ndim);
    for(int i= 0; i < n; i++)
        points.push_back(p);
    
    ZTable ztable;                      // genere une nouvelle table de permutation pour chaque realisation
    //~ initZTable(ztable, ndim/2);         // Zsampler genere des paires de dimensions
    randomZTable(ztable, ndim/2);         // Zsampler genere des paires de dimensions
    
    int i= 0;
    for(unsigned int y= 0; y < size; y++)
    for(unsigned int x= 0; x < size; x++, i++)
    {
        double tmp[2*nspp];
        for(int d= 0; d < ndim / 2; d++)
        {
            Z2D(ztable,
                x, y, // pixel
                1, nspp,                // spp
                Log2Int(size),          // 1 pixel dans une image size x size
                tmp,                    // spp samples pour la paire de dimensions d
                d);
            
            // re-ordonne les nspp samples
            for(int k= 0; k < nspp; k++)
            {
                points[i*nspp+k][2*d]= tmp[2*k];
                points[i*nspp+k][2*d+1]= tmp[2*k+1];
            }
        }
        
        //~ // force la stratification dans le pixel x, y sur une image size x size
        //~ points[i][0]= (x + points[i][0]) / double(size);
        //~ points[i][1]= (y + points[i][1]) / double(size);
    }
    
    assert(points.size() == n);
    return true;
}


const unsigned char digit_shuffle[24][4]= 
{
    {0, 1, 2, 3},
    {0, 1, 3, 2},
    {0, 2, 1, 3},
    {0, 2, 3, 1},
    {0, 3, 1, 2},
    {0, 3, 2, 1},
    {1, 0, 2, 3},
    {1, 0, 3, 2},
    {1, 2, 0, 3},
    {1, 2, 3, 0},
    {1, 3, 2, 0},
    {1, 3, 0, 2},
    {2, 0, 1, 3},
    {2, 0, 3, 1},
    {2, 1, 0, 3},
    {2, 1, 3, 0},
    {2, 3, 0, 1},
    {2, 3, 1, 0},
    {3, 0, 1, 2},
    {3, 0, 2, 1},
    {3, 1, 0, 2},
    {3, 1, 2, 0},
    {3, 2, 0, 1},
    {3, 2, 1, 0}
};

// permute l'arbre du code de morton
unsigned int shuffle_morton_tree( const unsigned int px, const unsigned int py, const int depth, const unsigned int seed= 0 )
{
    RNG rng;
    // chaque noeud re-ordonne ses fils
    rng.seed(seed);
    const unsigned char *shuffle= digit_shuffle[rng.sample_range(24)];
    
    unsigned int level= 1;
    unsigned int offset= 0;
    unsigned int code= 0;
    for(int i= depth -1; i >= 0; i--, level++)
    {
        unsigned int x= (px >> i) & 1;
        unsigned int y= (py >> i) & 1;
        unsigned int digit= y << 1 | x;
        
        // indice du noeud dans l'arbre de permutation
        // level 0 : 0                                                                      // base 0
        // level 1 : 00 | 01 | 02 | 03                                                      // base 1
        // level 2 : 000 001 002 003 | 010 011 012 013 | 020 021 022 023 | 030 031 032 033  // base 1+4= 5
        // level 3 :                                                                        // base 1+4+16= 21 = somme des puissances de 4
        // somme des puissances de 4 = \sum_{i=0}^{i<n} 4^i = (4^n -1)/3 = (2^{2n} -1)/3
        unsigned int base= ((1 << (2*level)) -1) / 3;
        assert(offset < (1<<(2*level)));
        unsigned int node= base + offset + digit;
        
        code= code << 2 | shuffle[digit];
        
        // position des fils
        offset= offset*4 + digit*4;
        
        // re-ordonne les fils
        rng.seed(seed ^ node);
        shuffle= digit_shuffle[rng.sample_range(24)];
    }
    
    return code;
}

// permute l'arbre du code de morton
unsigned int shuffle_morton_tree( const unsigned int px, const unsigned int py, const int depth, const int index, const int index_depth, const unsigned int seed= 0 )
{
    // etape 1 : permute les coordonnees du pixel
    // chaque noeud du quadtree / morton re-ordonne ses fils
    RNG rng;
    rng.seed(seed);
    const unsigned char *shuffle= digit_shuffle[rng.sample_range(24)];
    
    unsigned int level= 1;
    unsigned int offset= 0;
    unsigned int code= 0;
    for(int i= depth -1; i >= 0; i--, level++)
    {
        unsigned int x= (px >> i) & 1;
        unsigned int y= (py >> i) & 1;
        unsigned int digit= y << 1 | x;
        
        // indice du noeud dans l'arbre de permutation
        // level 0 : 0                                                                      // base 0
        // level 1 : 00 | 01 | 02 | 03                                                      // base 1
        // level 2 : 000 001 002 003 | 010 011 012 013 | 020 021 022 023 | 030 031 032 033  // base 1+4= 5
        // level 3 :                                                                        // base 1+4+16= 21 = somme des puissances de 4
        // somme des puissances de 4 = \sum_{i=0}^{i<n} 4^i = (4^n -1)/3 = (2^{2n} -1)/3
        unsigned int base= ((1u << (2*level)) -1) / 3;
        assert(offset < (1u<<(2*level)));
        unsigned int node= base + offset + digit;
        
        code= code << 2 | shuffle[digit];
        
        // position des fils
        offset= offset*4 + digit*4;
        
        // re-ordonne les fils
        rng.seed(seed ^ node);
        shuffle= digit_shuffle[rng.sample_range(24)];
    }
    
    // etape 2 : permute aussi l'index du sample
    assert((index_depth & 1) == 0);
    
    for(int i= index_depth -2; i >= 0; i-=2, level++)
    {
        unsigned int digit= (index >> i) & 3;
        
        unsigned int base= ((1u << (2*level)) -1) / 3;
        assert(offset < (1u<<(2*level)));
        unsigned int node= base + offset + digit;
        
        code= code << 2 | shuffle[digit];
        
        // position des fils
        offset= offset*4 + digit*4;
        
        // re-ordonne les fils
        rng.seed(seed ^ node);
        shuffle= digit_shuffle[rng.sample_range(24)];
    }
    
    return code;
}


// hierarchical morton shuffle, cf Z sampler
bool sample_morton( std::vector< VecXDynamic >& points, const int n, const int ndim, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    // verifie que spp est un carre 
    int size= std::sqrt(n);
    if(size*size != n)
        return false;   // pas possible
    
    VecXDynamic sample(ndim);
    
    std::random_device hwseed;
    unsigned int morton_seed= hwseed();
    int morton_depth= std::log2(size);
    
    for(int y= 0; y < size; y++)
    for(int x= 0; x < size; x++)
    {
        unsigned int index= shuffle_morton_tree(x, y, morton_depth, morton_seed);
        for(int d= 0; d < ndim; d++)
            sample[d] = double(sobols[d].getSobolInt(index)) / double(UINT32SOBOLNORM);
        
        //~ // force la stratification dans le pixel x, y sur une image size x size
        //~ sample[0]= (x + sample[0]) / double(size);
        //~ sample[1]= (y + sample[1]) / double(size);
        
        points.push_back( sample );
    }
    
    assert(points.size() == n);
    return true;
}

// hierarchical morton shuffle, cf Z sampler,
// n'utilise que les dimensions 01 + permutation morton differente par paire
bool sample_morton01( std::vector< VecXDynamic >& points, const int n, const int ndim, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    // verifie que spp est un carre 
    int size= std::sqrt(n);
    if(size*size != n)
        return false;   // pas possible
    
    // une permutation par paire de dimensions
    std::random_device hwseed;
    std::vector<unsigned int> morton_seeds;
    for(int pair= 0; pair < ndim / 2; pair++)
        morton_seeds.push_back( hwseed() );

    VecXDynamic sample(ndim);
    
#if 0
    // permute les coordonnees du pixel
    int morton_depth= std::log2(size);
    for(int y= 0; y < size; y++)
    for(int x= 0; x < size; x++)
    {
        for(int pair= 0; pair < ndim / 2; pair++)
        {
            //~ unsigned int index= shuffle_morton_tree(x, y, morton_depth, morton_seeds[pair]);
            // sobol shift...
            unsigned int index= shuffle_morton_tree(x, y, morton_depth, morton_seeds[pair]) + n;
            
            sample[2*pair] = double(sobols[0].getSobolInt(index)) / double(UINT32SOBOLNORM);
            sample[2*pair+1] = double(sobols[1].getSobolInt(index)) / double(UINT32SOBOLNORM);
        }
        
        //~ // force la stratification dans le pixel x, y sur une image size x size
        //~ sample[0]= (x + sample[0]) / double(size);
        //~ sample[1]= (y + sample[1]) / double(size);
        
        points.push_back( sample );
    }
#else
    // permute les coordonnees du pixel (0, 0) et de l'indice du sample
    int index_depth= std::log2(n);
    for(int i= 0; i < n; i++)
    {
        for(int pair= 0; pair < ndim / 2; pair++)
        {
            unsigned int index= shuffle_morton_tree(0, 0, 1, i, index_depth, morton_seeds[pair]);
            
            sample[2*pair] = double(sobols[0].getSobolInt(index)) / double(UINT32SOBOLNORM);
            sample[2*pair+1] = double(sobols[1].getSobolInt(index)) / double(UINT32SOBOLNORM);
        }
        
        points.push_back( sample );
    }
#endif

    assert(points.size() == n);
    return true;
}

// hierarchical morton shuffle, cf Z sampler
bool sample_morton( std::vector< VecXDynamic >& points, const int size, const int n, const int ndim, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    //~ int size= 4;
    int nspp= n / (size * size);
    if(nspp != RoundUpPow2(nspp))
        return false;
    if(nspp*size*size != n)
        return false;
    
    VecXDynamic sample(ndim);
    
    std::random_device hwseed;
    unsigned int morton_seed= hwseed();
    int morton_depth= std::log2(size);
    
    for(int y= 0; y < size; y++)
    for(int x= 0; x < size; x++)
    {
        unsigned int index= shuffle_morton_tree(x, y, morton_depth, morton_seed);
        for(int i= 0; i < nspp; i++)
        {
            for(int d= 0; d < ndim; d++)
                sample[d] = double(sobols[d].getSobolInt(index*nspp + i)) / double(UINT32SOBOLNORM);
            
            //~ // force la stratification dans le pixel x, y sur une image size x size
            //~ sample[0]= (x + sample[0]) / double(size);
            //~ sample[1]= (y + sample[1]) / double(size);
            
            points.push_back( sample );
        }
    }
    
    assert(points.size() == n);
    return true;
}


// random shuffle
bool sample_morton_white( std::vector< VecXDynamic >& points, const int n, const int ndim, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    // verifie que n est un carre 
    int size= std::sqrt(n);
    if(size*size != n)
        return false;   // pas possible
    
    VecXDynamic sample(ndim);
    
    std::random_device hwseed;
    unsigned int morton_seed= hwseed();
    int morton_depth= std::log2(size);
    
    std::vector<uint32_t> sobol_indices;
    for(int i= 0; i < size*size; i++)
        sobol_indices.push_back(i);
    
    std::default_random_engine shuffle(hwseed());
    std::shuffle(sobol_indices.begin(), sobol_indices.end(), shuffle);
    
    for(int y= 0; y < size; y++)
    for(int x= 0; x < size; x++)
    {
        unsigned int index= sobol_indices[y*size +x];
        for(int d= 0; d < ndim; d++)
            sample[d] = double(sobols[d].getSobolInt(index)) / double(UINT32SOBOLNORM);
        
        //~ // force la stratification dans le pixel x, y sur une image size x size
        //~ sample[0]= (x + sample[0]) / double(size);
        //~ sample[1]= (y + sample[1]) / double(size);
        
        points.push_back( sample );
    }
    
    assert(points.size() == n);
    return true;
}

// random shuffle
bool sample_morton_white( std::vector< VecXDynamic >& points, const int size, const int n, const int ndim, const std::vector< SobolGenerator1D<uint32_t> >& sobols )
{
    assert(!sobols.empty());
    points.clear();
    int nspp= n / (size * size);
    if(nspp != RoundUpPow2(nspp))
        return false;
    if(nspp*size*size != n)
        return false;
    
    VecXDynamic sample(ndim);
    
    std::random_device hwseed;
    unsigned int morton_seed= hwseed();
    int morton_depth= std::log2(size);
    
    std::vector<uint32_t> sobol_indices;
    for(int i= 0; i < size*size; i++)
        sobol_indices.push_back(i);
    
    std::default_random_engine shuffle(hwseed());
    std::shuffle(sobol_indices.begin(), sobol_indices.end(), shuffle);
    
    for(int y= 0; y < size; y++)
    for(int x= 0; x < size; x++)
    {
        unsigned int index= sobol_indices[y*size +x];
        for(int i= 0; i < nspp; i++)
        {
            for(int d= 0; d < ndim; d++)
                sample[d] = double(sobols[d].getSobolInt(index*nspp + i)) / double(UINT32SOBOLNORM);
            
            //~ // force la stratification dans le pixel x, y sur une image size x size
            //~ sample[0]= (x + sample[0]) / double(size);
            //~ sample[1]= (y + sample[1]) / double(size);
            
            points.push_back( sample );
        }
    }
    
    assert(points.size() == n);
    return true;
}


// pbrt (0,2) seq sampler
extern
void ZTSequence2D(
    int nSamplesPerPixelSample,
    int nPixelSamples,
    std::vector< std::vector<Point2f> >& samples2D,
    int nDimensions,
    unsigned int seed    
);

bool sample_zerotwo( std::vector< VecXDynamic >& points, const int spp, const int ndim )
{
    // prealloue tous les samples
    points.clear();
    VecXDynamic p(ndim);
    for(int i= 0; i < spp; i++)
        points.push_back(p);

    // temporaire
    std::vector< std::vector<Point2f> > tmp;
    // pre-alloue les samples ZT
    for(int d= 0; d < ndim/2; d++)
        tmp.push_back( std::vector<Point2f>(spp) );

    std::random_device hwseed;
    unsigned int seed= hwseed();
    
    ZTSequence2D(1, spp, tmp, tmp.size(), seed);

    // re-organise les points
    for(int i= 0; i < spp; i++)
    for(int d= 0; d < ndim/2; d++)
    {
        points[i][2*d]= tmp[d][i].x;
        points[i][2*d+1]= tmp[d][i].y;
    }

    return true;
}


double integrate( const std::vector< VecXDynamic >& points )
{
    assert(!points.empty());
#ifdef GAUSS
    return calculate_mse(points, 2, 16384); // Gauss
#else
    return calculate_mse(points, 1, 16384); // Heaviside
#endif
}

double discrepancy( const std::vector< VecXDynamic >& points )
{
    assert(!points.empty());
    return generalizedL2(points, points[0].dim());
}


double integrate( const std::vector<int>& integrands, const std::vector< VecXDynamic >& points )
{
    assert(points.size() > 0);
    assert(points[0].dim() == 4);
    
//~ #ifdef GAUSS
    //~ #error can't use GAUSS and rank integrands...
//~ #endif
    
    double mse_sum= 0;
    const int n= int(integrands.size());
    for(int i= 0; i < n; i++)
    {
        int k= integrands[i];
        
    //~ #ifdef GAUSS
        //~ const double integral= tab_Gauss4D[k].integral;
        //~ const double *mu= tab_Gauss4D[k].mu;
        //~ const double *mxCInv= tab_Gauss4D[k].mxCInv;
    //~ #else
        const double integral= tab_Heaviside4D[k].integral;
        const double *mu= tab_Heaviside4D[k].muDiscotinuity;
        const double *normal= tab_Heaviside4D[k].normal;
    //~ #endif
        
        // estimate 
        double sum= 0;
        for(int p= 0; p < int(points.size()); p++)
        {
        //~ #ifdef GAUSS
            //~ sum+= getMultivariateGaussian(4, points[p].data(), mu, mxCInv);
            
        //~ #else
            double dot= 0;
            for(int d= 0; d < 4; d++)
                dot+= (points[p][d] - mu[d]) * normal[d];
            
            sum+= (dot > 0) ? 1 : 0;
        //~ #endif
        }
        sum/= double(points.size());
        
        // mse
        double mse= (integral - sum) * (integral - sum);
        mse_sum+= mse;
    }
    mse_sum/= n;
    return mse_sum;
}

#ifdef RANK
// moyenne des rangs des fonctions de tests
double rank( const std::vector< VecXDynamic >& points )
{
    assert(points.size() > 0);
    assert(points[0].dim() == 4);
    int n= int(points.size());
    
    std::vector<int> ranks(16384);
    
    double mranks= 0;
    #pragma omp parallel for reduction(+:mranks)
    for(int k= 0; k < 16384; k++)
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix(n, n);
        matrix.Zero(n, n);
        
    #ifdef GAUSS
        // gauss
        const double *mu = tab_Gauss4D[k].mu;
        const double *mxCInv = tab_Gauss4D[k].mxCInv;
        //~ const double mu[4]= { 0.5, 0.5, 0.5, 0.5 };
        //~ const double mxCInv[16]= { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
        
        for(int i= 0; i < n; i++)
        {
            VecX<4> point;
            point[0]= points[i][0];
            point[1]= points[i][1];
            
            for(int j= 0; j < n; j++)
            {
                point[2]= points[j][2];
                point[3]= points[j][3];
                
                // eval function
                matrix(i, j)= getMultivariateGaussian(4, point.data(), mu, mxCInv);
            }
        }
        
    #else
        // heaviside
        const double *mu= tab_Heaviside4D[k].muDiscotinuity;
        const double *normal= tab_Heaviside4D[k].normal;
        
        for(int i= 0; i < n; i++)
        {
            VecX<4> point;
            point[0]= points[i][0];
            point[1]= points[i][1];
            
            for(int j= 0; j < n; j++)
            {
                point[2]= points[j][2];
                point[3]= points[j][3];
                
                // eval function
                double dot= 0;
                for(int d= 0; d < 4; d++)
                    dot+= (point[d] - mu[d]) * normal[d];
                
                matrix(i, j)= (dot > 0) ? 1 : 0;
            }
        }
    #endif
        
        // compute rank
        Eigen::FullPivLU< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > LU(matrix);
        LU.setThreshold(1e-5);
        int rank= LU.rank();
        
        ranks[k]= rank;
        mranks+= double(rank) / double(16384);
    }
    
    // export ranks
    {
        std::vector<int> sorted= ranks;
        std::sort(sorted.begin(), sorted.end());
        
    #ifdef GAUSS
        std::ofstream ofs("ranks_gauss.dat", std::ios_base::out);
    #else
        std::ofstream ofs("ranks_heaviside.dat", std::ios_base::out);
    #endif
        
        assert(bool(ofs));
        
        for(int i= 0; i < int(ranks.size()); i++)
        {
            ofs << ranks[i] << " " << sorted[i] << std::endl;
        }
        
        ofs.close();
    }
    
    // export histogram
    {
        std::vector<double> histogram(256, 0);
        for(int i= 0; i < int(ranks.size()); i++)
        {
            int rank= int(ranks[i]);
            if(rank >= int(histogram.size()))
                histogram.resize(rank+1);
            
            histogram[rank]+= 1.0 / double(ranks.size());
        }
        
        double p= 0;
        for(int i= 0; i < int(histogram.size()); i++)
            p+= histogram[i];
        std::cout << "histogram " << int(histogram.size()) << " bins, normalization " << p << std::endl;
        
    #ifdef GAUSS
        std::ofstream ofs("histogram_gauss.dat", std::ios_base::out);
    #else
        std::ofstream ofs("histogram_heaviside.dat", std::ios_base::out);
    #endif
        assert(bool(ofs));
        
        for(int i= 0; i < int(histogram.size()); i++)
        {
            ofs << histogram[i] << std::endl;
        }
        
        ofs.close();
    }
    
    return mranks;
}
#endif

#ifndef RANK
// selectionne une distribution de fonctions heaviside 4d en fonction d'un histogramme de rangs
std::vector<int> match_integrands( const std::string& histogram_filename )
{
    std::vector<double> histogram;
    {
        std::cout << "reading '" << histogram_filename << "'...\n";
        std::ifstream ifs(histogram_filename, std::ios_base::in);
        if(!ifs.is_open())
        {
            std::cout << "[error] reading rank histogram...\n";
            exit(1);
        }
        
        double p;
        while(ifs >> p)
            histogram.push_back(p);
        ifs.close();
        
        std::cout << "histogram: " << histogram.size() << std::endl;
    }
    
    std::vector<int> ranks;
    std::vector<double> ranks_histogram(256);
    std::vector<double> ranks_cdf;
    {
        std::cout << "reading 'ranks_heaviside.dat'...\n";
        std::ifstream ifs("ranks_heaviside.dat", std::ios_base::in);
        if(!ifs.is_open())
        {
            std::cout << "[error] reading ranks...\n";
            exit(1);
        }
        
        int rank, sorted;
        while(ifs >> rank >> sorted)
            ranks.push_back(rank);
        ifs.close();
        
        std::cout << "ranks: " << ranks.size() << std::endl;
        
        for(int i= 0; i < int(ranks.size()); i++)
        {
            int rank= int(ranks[i]);
            if(rank >= int(ranks_histogram.size()))
                ranks_histogram.resize(rank+1);
            
            ranks_histogram[rank]+= 1.0 / double(ranks.size());
        }
        
        // match distributions
        double sum= 0;
        for(int i= 0; i < int(ranks.size()); i++)
        {
            int rank= ranks[i];
            double p= histogram[rank] / ranks_histogram[rank];
            sum+= p;
            ranks_cdf.push_back(sum);
        }
        
        for(int i= 0; i < ranks_cdf.size(); i++)
            ranks_cdf[i]/= sum;
    }
    
    // resample
    const unsigned int seed= 1;
    RNG rng(seed);
    
    std::vector<int> integrands;
    for(int i= 0; i < 16384; i++)
    {
        double x= rng.sample_double();
        // invert cdf
        int k= 0;
        for(; k < int(ranks_cdf.size()); k++)
            if(x < ranks_cdf[k])
                break;
        
        integrands.push_back(k);
    }
    
    // export distribution 
    {
        std::vector<double> histogram(256, 0);
        for(int i= 0; i < int(integrands.size()); i++)
        {
            int k= integrands[i];
            int rank= int(ranks[k]);
            if(rank >= int(histogram.size()))
                histogram.resize(rank+1);
	
            histogram[rank]+= 1.0 / double(integrands.size());
        }
        
        std::ofstream ofs("distribution.dat", std::ios_base::out);
        assert(bool(ofs));
        
        for(int i= 0; i < int(histogram.size()); i++)
        {
            ofs << histogram[i] << std::endl;
        }
        
        ofs.close();
    }
    
    return integrands;
}
#endif

struct lessVecX
{
    bool operator() ( const VecXDynamic& a, const VecXDynamic& b ) const
    {
        assert(a.dim() == b.dim());
        const int n= a.dim();
        for(int d= 0; d < n; d++)
            if(a[d] != b[d])
                return a[d] < b[d];
        
        return false;
    }
};

int main( int argc, char **argv )
{
    CLI::App app { "Integration test in nD using sobol++."};

    int ndim= 2;
    app.add_option("-d", ndim, "Dimension (defaukt 2)");
    std::string dimensions_string;
    app.add_option("--dims", dimensions_string, "explicit dimensions \"i j k \"");
    
    std::string dir_vectors_fname= "";
    app.add_option("--dirs", dir_vectors_fname, "File name of the Sobol intialization table (e.g. ../../../data/sobol_init_tab.dat)");
    std::string output_fname;
    app.add_option("-o,--output", output_fname, "Output filename for the stats (ascii .Dat).")->required();
    int realizations= 1;
    app.add_option("-n,--nbReal", realizations, "Number of realizations (with random owen++scrambling) (def 1)");
    int maxK= 16;
    app.add_option("-m,--maxK", maxK, "K for the maximum number of spp (2^K, def=16)");
    
    int owen_depth= 32;
    //~ app.add_option("--owen_depth", owen_depth, "Owen++ Permutation depth, min(32, log2(N)+depth), default 32");
    app.add_option("--owen_depth", owen_depth, "Owen++ Permutation depth, default 32");
    
    bool disableDiscrepancy= false;
    //~ bool disableDiscrepancy= true;
    app.add_flag("--disableDiscrepancy", disableDiscrepancy, "Disable the discrepancy test");

    std::string match_integrands_filename;
    app.add_option("--ranks", match_integrands_filename, "selects Heaviside functions, requires '-d 4' and #undef GAUSS");
    
    std::string export_prefix;
    app.add_option("--export", export_prefix, "Output prefix for the points 'prefixXXXX.dat'.");

    bool export_sort= false;
    app.add_flag("--sort", export_sort, "sort points, requires --export");
    
    bool sobol= false;
    app.add_flag("--sobol", sobol, "Sobol test");
    bool sobol_offset= false;
    app.add_flag("--offset", sobol_offset, "Sobol test");
    bool csobol= false;
    app.add_flag("--csobol", csobol, "Cascaded Sobol test");
    bool owen= false;
    app.add_flag("--owen", owen, "Owen test");
    bool owenpp= false;
    app.add_flag("--owenpp", owenpp, "Owen++ test");
    bool owenpp01= false;
    app.add_flag("--owenpp01", owenpp01, "Owen++ 01 test");
    bool zsampler= false;
    app.add_flag("--zsampler", zsampler, "ZSampler test");
    bool morton= false;
    app.add_flag("--morton", morton, "Morton shuffle test");
    bool morton01= false;
    app.add_flag("--morton01", morton01, "Morton 01 shuffle test");
    bool white= false;
    app.add_flag("--white", white, "Morton White shuffle test");
    bool zerotwo= false;
    app.add_flag("--zerotwo", zerotwo, "Pbrt (0,2) seq test");
    bool fromfile= false;
    app.add_flag("--fromfile", zerotwo, "Import samples from files");
    std::string prefix;
    app.add_flag("--prefix", prefix, "Filename prefix for point set import, filename = prefix${nbpts}.dat");
    int zeropadding=0;
    app.add_flag("--padding", zeropadding, "Size of zero padding for sample filenames");
  
    int size= 0;
    app.add_option("--size", size, "image size");
  
    CLI11_PARSE(app, argc, argv)

    std::vector< SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    if(!dir_vectors_fname.empty())
    {
        loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file and fill appropriate structures
        
        // debug
	    for(int i= 0; i < ndim; i++) 
        {
	        std::cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
	        for(int j = 0; j < sobols[i].s; j++) 
                std::cout << " " << sobols[i].m[j];
	        std::cout  << std::endl;
	    }        
    }

    std::vector<int> dimensions;
    for(int i= 0; i < ndim; i++)
        dimensions.push_back(i);

    if(!dimensions_string.empty())
    {
        // parse explicit dimensions
        std::istringstream in(dimensions_string);
        for(int i= 0; i < ndim; i++)
            if(!(in >> dimensions[i]))
                break;
        
        std::cout << "using explicit dimensions: ";
        for(int i= 0; i < ndim; i++)
            std::cout << dimensions[i] << " ";
        std::cout << std::endl;
    }

    std::vector<int> integrands;
#ifndef RANK
    if(!match_integrands_filename.empty())
    {
        //~ integrands= match_integrands("cornell_diffuse_ranks.dat");
        integrands= match_integrands(match_integrands_filename);
        ndim= 4;
        std::cout << "using 4d rank integrands...\n";
    }
    else
#endif
    {
    #ifdef GAUSS
        std::cout << "using gauss functions...\n";
    #else
        std::cout << "using heaviside functions...\n";
    #endif
    }
    
    std::ofstream ofs(output_fname, std::ios_base::out);
    assert(bool(ofs));
    
    for(int k= 1; k <= maxK; k++)
    {
        size_t spp= std::pow(2, k);
        std::cout << "### " << spp << std::endl;
        
        std::vector< VecXDynamic > points;
        
        int results= 0;
        double integration_result= 0;
        double discrepancy_result= 0;
        double rank_result= 0;

        std::string filename = std::to_string(zeropadding);
        filename = std::string( std::max(0,zeropadding -  int(filename.length())), '0').append(filename);
        std::ifstream in (filename);

        for(int n= 0; n < realizations; n++)
        {
            bool samples= false;
            
            if(sobol) 
                //~ samples= sample_sobol(points, spp, ndim, sobols);
                samples= sample_sobol(points, spp, dimensions, sobols);
            else if(sobol_offset) 
                samples= sample_sobol_offset(points, spp, ndim, sobols);
            else if(csobol)
                samples= sample_cascaded_sobol(points, spp, ndim, sobols, owen_depth);
            else if(owen)
                samples= sample_owen(points, spp, dimensions, sobols);
            else if(owenpp)
                samples= sample_owenpp(points, spp, dimensions, sobols, owen_depth);
            else if(owenpp01)
                samples= sample_owenpp01(points, spp, ndim, sobols);
            
            else if(zsampler && size == 0)
                samples= sample_zsampler(points, spp, ndim);            // 1 sample per pixel
            else if(zsampler && size != 0)
                samples= sample_zsampler(points, size, spp, ndim);      // size x size image, samples per pixel
            
            else if(morton && size == 0)
                samples= sample_morton(points, spp, ndim, sobols);
            else if(morton && size != 0)
                samples= sample_morton(points, size, spp, ndim, sobols);
            else if(morton01 && size == 0)
                samples= sample_morton01(points, spp, ndim, sobols);
            
            else if(white && size == 0)
                samples= sample_morton_white(points, spp, ndim, sobols);
            else if(white && size != 0)
                samples= sample_morton_white(points, size, spp, ndim, sobols);
            
            else if(zerotwo)
                samples= sample_zerotwo(points, spp, ndim);
            else if(fromfile){
                samples = not in.eof();
                points.clear();
                for (int idPt = 0; idPt < n; ++idPt){
                    points.emplace_back(ndim);
                    for (int idComp = 0; idComp < ndim; ++ndim){
                        in >> points.back()[idComp];
                    }
                }
                char trash;
                in >> trash;
            }
            
            else
            {
                printf("[error] no sampler...\n");
                break;
            }
            
            if(samples)
            {
                results++;
                assert(!points.empty());
                
            #ifndef RANK
                if(!integrands.empty())
                    // n'utilise que les fonctions selectionnees
                    integration_result+= integrate(integrands, points);
                else
                    integration_result+= integrate(points);
                
                if(!disableDiscrepancy)
                    discrepancy_result+= discrepancy(points);
            #endif
            
            #ifdef RANK
                rank_result+= rank(points);
            #endif
                
                if(!export_prefix.empty() && n == 0)
                {
                    if(export_sort)
                        std::sort(points.begin(), points.end(), lessVecX());
                    
                    std::string fname= export_prefix + std::to_string(spp) + ".dat";
                    std::ofstream ofs(fname, std::ios_base::out);
                    for(int i= 0; i < spp; i++)
                    {
                        for(int d=0; d < ndim; d++)
                            ofs << points[i][d] << " ";
                        ofs << std::endl;
                    }
                    ofs.close();
                }
            }
        }
        
        if(results)
        {
            integration_result/= double(results);
            discrepancy_result/= double(results);
            
        #ifdef RANK
            rank_result/= double(results);
            std::cout << "rank " << rank_result << std::endl;
        #endif
            
            ofs << spp << " " << integration_result << " ";
            if(!disableDiscrepancy)
                ofs << discrepancy_result;
            
            ofs << std::endl;
        }
    }
    ofs.close();
    
    return 0;
}

