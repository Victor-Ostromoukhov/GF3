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

#include <CLI11.hpp>

#include "../Points/VecX.h"
#include "../Samplers/Random.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

#include "../Screenspace/ZSampler/z.h"


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

int main( int argc, char **argv )
{
    CLI::App app { "output points to dat files"};

    std::string dimensions_string;
    app.add_option("--dims", dimensions_string, "explicit dimensions \"i j k ...\"");
    int ndim= 4;
    app.add_option("-d", ndim, "Dimension (default 4)");
    std::string dir_vectors_fname= "";
    app.add_option("--dirs", dir_vectors_fname, "File name of the Sobol intialization table (e.g. ../../../data/sobol_init_tab.dat)");
    std::string output_fname;
    app.add_option("-o,--output", output_fname, "Output filename for the points (ascii .Dat).")->required();
    int realizations= 1;
    app.add_option("-m,--nbReal", realizations, "Number of realizations (with random owen++scrambling) (def 1)");
    int n = 4;
    app.add_option("-n,--nbPts", n, "Number of Points");

    bool sobol= false;
    app.add_flag("--sobol", sobol, "Sobol test");
    bool sobol_offset= false;
    app.add_flag("--offset", sobol_offset, "Sobol test");
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

    CLI11_PARSE(app, argc, argv)

    std::vector< SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    if(!dir_vectors_fname.empty())
        loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file and fill appropriate structures

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

    std::ofstream ofs(output_fname, std::ios_base::out);
    assert(bool(ofs));
    for(int idReal = 0; idReal < realizations; idReal++) {
            bool samples= false;

            std::vector<VecXDynamic> points;

            if(sobol)
                //~ samples= sample_sobol(points, n, ndim, sobols);
                samples= sample_sobol(points, n, dimensions, sobols);
            else if(sobol_offset)
                samples= sample_sobol_offset(points, n, ndim, sobols);
            else if(owen)
                samples= sample_owen(points, n, dimensions, sobols);
            else if(owenpp)
                samples= sample_owenpp(points, n, dimensions, sobols, 32);
            else if(owenpp01)
                samples= sample_owenpp01(points, n, ndim, sobols);

            else if(zsampler)
                samples= sample_zsampler(points, n, ndim);            // 1 sample per pixel

            else if(morton)
                samples= sample_morton(points, n, ndim, sobols);
            else if(morton01)
                samples= sample_morton01(points, n, ndim, sobols);

            else if(white)
                samples= sample_morton_white(points, n, ndim, sobols);

            else if(zerotwo)
                samples= sample_zerotwo(points, n, ndim);


            else
            {
                printf("[error] no sampler...\n");
                break;
            }

            if(samples)
            {
                if (idReal != 0){
                    ofs << '#' << std::endl;
                }
                for (auto& v : points){
                    ofs << v << std::endl;
                }
            }
    }
    ofs.close();

    return 0;
}

