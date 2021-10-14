
#ifndef _SOBOLPP_SAMPLER_H
#define _SOBOLPP_SAMPLER_H 

#include "real.h"
#include "sampler.h"

#include "../includes/SeedTileHelper.hpp"
#include "../Samplers/BlueTile.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"
#include "../Samplers/SobolCascade.h"


struct SobolPPSampler : public Sampler
{
    SobolPPSampler( ) : Sampler(), 
        m_shared_seeds(nullptr), m_shared_sobols(nullptr), m_dimensions(), m_owen_bits(0), m_tile_size(0), m_width(0), m_height(0), m_spp(0) 
    {}
    
    explicit SobolPPSampler( const std::vector<int>& d, const int width, const int height, const size_t spp, 
        const char *seed_filename, const char *sobol_filename, const int owen_bits, const unsigned int s ) : Sampler(d.size(), width*height*spp, s), 
        m_shared_seeds(nullptr), m_shared_sobols(nullptr), m_dimensions(d), m_owen_bits(32),
        m_tile_size(0), m_width(width), m_height(height), m_spp(spp)
    {
        if(seed_filename == nullptr)
        {
            // init random, sans chargement de fichier
            int tile_size= std::min(m_width, m_height);
            //~ int tile_size= 1;   // sequence coherente / identique par pixel
            printf("[sobolpp] random init %dx%d %dd\n", tile_size, tile_size, m_ndim);
            
            auto *seeds= new BlueTile(tile_size, m_ndim, 1);
            seeds->random_init(s);
            
            m_shared_seeds= seeds;
        }
        else
        {
            auto *seeds= new BlueTile();
            seeds->read(seed_filename);
            assert(seeds->size());
            
            m_shared_seeds= seeds;
        }
        
        m_tile_size= m_shared_seeds->width();
        if(m_shared_seeds->dimensions() != m_ndim)
        {
            printf("[error] seed tile %dd != sampler %dd\n", m_shared_seeds->dimensions(), m_ndim);
            exit(1);
        }
        
        auto *sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *sobols);    // read sobols from file and fill appropriate structures
        
        m_shared_sobols= sobols;
        int sobol_bits= std::log2(m_spp);
        m_owen_bits= std::min(sobol_bits + owen_bits, 32);
    }
    
    SobolPPSampler( const SobolPPSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_shared_seeds(sampler.m_shared_seeds), m_shared_sobols(sampler.m_shared_sobols), m_dimensions(sampler.m_dimensions), m_owen_bits(sampler.m_owen_bits),
        m_tile_size(sampler.m_tile_size), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp) 
    {}
    // tous les threads, referencent la meme copie des seeds / des sobols.
    // utiliser sampler::release() pour la detruire, uniquement sur main_sampler, pas sur les samplers qui ne detiennent que la reference
    
    SobolPPSampler& operator= ( const SobolPPSampler& ) = delete;
    SobolPPSampler( const SobolPPSampler& ) = delete;
    
    ~SobolPPSampler( ) 
    {
        printf("sobolpp max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("sobolpp release...\n");
        delete m_shared_seeds;
        delete m_shared_sobols;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new SobolPPSampler(*this, s);
    }
    
    void seed( const unsigned int s ) 
    { 
        // pas de seeds aleatoires 
    }
    
    void index( const size_t id )
    {
        m_maxdim= std::max(m_maxdim, m_dim);
        
        m_dim= 0;
        m_index= id;
        
        // associe un pixel a l'indice
        int pi= id / m_spp;
        int px= pi % m_width;
        int py= pi / m_width;
        
        assert(px >= 0 && px < m_width);
        assert(py >= 0 && py < m_height);
        
        // recupere le sample
        assert(m_data.size());
        assert(m_data.size() == m_dimensions.size());
        // uniquement les dimensions necessaires...
        for(int d= 0; d < m_data_ndim; d++)
        {
            assert(m_dimensions[d] != -1);
            unsigned int owen_seed= (*m_shared_seeds)(px % m_tile_size, py % m_tile_size, d, 0);
            m_data[d]= getOwenPlus1D_with_seed(*m_shared_sobols, m_dimensions[d], id % m_spp, owen_seed, m_owen_bits);
        }
        
        // puis remappe les dimensions...
        m_values[0]= Float(0.5);    // centre du pixel
        m_values[1]= Float(0.5);
        m_values[2]= Float(0.5);    // centre de la lentille
        m_values[3]= Float(0.5);
        for(int i= 4; i < m_ndim; i++) 
            m_values[i]= Float(0);   // force un zero dans les dimensions non utilisees
            
        for(int i= 0; i < m_ndim; i++)
            if(m_data_remap[i] != -1)
                m_values[i]= m_data[m_data_remap[i]];
        
        // force une stratification sur les dimensions xy
        m_values[0]= (Float(px) + m_values[0]) / Float(m_width);
        m_values[1]= (Float(py) + m_values[1]) / Float(m_height);        
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);      // ne doit pas arriver, c'est index() qui fait tout le boulot...
    }

protected:
    const BlueTile *m_shared_seeds;
    const std::vector< SobolGenerator1D<uint32_t> > *m_shared_sobols;
    std::vector<int> m_dimensions;
    
    int m_owen_bits;
    int m_tile_size;
    int m_width;
    int m_height;
    int m_spp;
};


struct CSobolPPSampler : public Sampler
{
    CSobolPPSampler( ) : Sampler(), 
        m_shared_seeds(nullptr), m_shared_sobols(nullptr),
        m_owen_bits(0), m_sobol_bits(0), 
        m_tile_size(0), m_width(0), m_height(0), m_spp(0) 
    {}
    
    explicit CSobolPPSampler( const int ndim, const int width, const int height, const size_t spp, 
        const char *seed_filename, const char *sobol_filename, const int owen_bits, const unsigned int s ) : Sampler(ndim, width*height*spp, s), 
        m_shared_seeds(nullptr), m_shared_sobols(nullptr), 
        m_owen_bits(32), m_sobol_bits(0),
        m_tile_size(0), m_width(width), m_height(height), m_spp(spp)
    {
        if(seed_filename == nullptr)
        {
            // init random, sans chargement de fichier
            int tile_size= std::min(m_width, m_height);
            //~ int tile_size= 1;   // sequence coherente / identique par pixel
            printf("[csobol] random init %dx%d %dd\n", tile_size, tile_size, m_ndim);
            
            auto *seeds= new BlueTile(tile_size, m_ndim, 1);
            seeds->random_init(s);
            
            m_shared_seeds= seeds;
            
            //~ m_shared_seeds->write("seeds_export.dat");
        }
        else
        {
            auto *seeds= new BlueTile();
            seeds->read(seed_filename);
            assert(seeds->size());
            
            m_shared_seeds= seeds;
        }
        
        m_tile_size= m_shared_seeds->width();
        if(m_shared_seeds->dimensions() != m_ndim)
        {
            printf("[error] seed tile %dd != sampler %dd\n", m_shared_seeds->dimensions(), m_ndim);
            exit(1);
        }
        
        auto *sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *sobols);    // read sobols from file and fill appropriate structures
        
        m_shared_sobols= sobols;
        m_sobol_bits= std::log2(m_spp);
        m_owen_bits= std::min(m_sobol_bits + owen_bits, 32);
        
        if(1)
        {
            for(int i= 0; i < m_ndim; i++) 
            {
                printf("%d : %u %u %u ", i, (*m_shared_sobols)[i].d, (*m_shared_sobols)[i].s, (*m_shared_sobols)[i].a);
                for(int j = 0; j < (*m_shared_sobols)[i].s; j++)
                    printf(" %u", (*m_shared_sobols)[i].m[j]);
                printf("\n");
            }
            
            printf("[csobol] ndim %dd, tile %d, spp %d, sobol bits %d, owen bits %d\n", 
                m_ndim, m_tile_size, m_spp,
                m_sobol_bits, m_owen_bits);
            
            if(0)
            {
                char tmp[1024];
                sprintf(tmp, "export_cascade_%05d.dat", m_spp);
                FILE *out= fopen(tmp, "wt");
                
                std::vector<Float> sample(m_ndim);
                for(int y= 0; y < 8; y++)
                for(int x= 0; x < 8; x++)
                {
                    for(int i= 0; i < m_spp; i++)
                    {
                        const uint32_t *seeds= (*m_shared_seeds)(x, y, 0);
                        // genere les n dimensions du sample
                        getUniformND(sample.data(), *m_shared_sobols, seeds, m_ndim, i, m_sobol_bits, m_owen_bits, true);    // owen sur log2(n) + 6
                        
                        for(int d= 0; d < m_ndim; d++)
                            fprintf(out, "%.10f ", sample[d]);
                        fprintf(out, "\n");
                    }
                    
                    fprintf(out, "#\n");
                }
                
                fclose(out);
            }
        }
    }
    
    CSobolPPSampler( const CSobolPPSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_shared_seeds(sampler.m_shared_seeds), m_shared_sobols(sampler.m_shared_sobols), m_owen_bits(sampler.m_owen_bits), m_sobol_bits(sampler.m_sobol_bits),
        m_tile_size(sampler.m_tile_size), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp) 
    {}
    // tous les threads, referencent la meme copie des seeds / des sobols.
    // utiliser sampler::release() pour la detruire, uniquement sur main_sampler, pas sur les samplers qui ne detiennent que la reference
    
    CSobolPPSampler& operator= ( const CSobolPPSampler& ) = delete;
    CSobolPPSampler( const CSobolPPSampler& ) = delete;
    
    ~CSobolPPSampler( ) 
    {
        printf("cascaded sobol max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("cascaded sobol release...\n");
        delete m_shared_seeds;
        delete m_shared_sobols;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new CSobolPPSampler(*this, s);
    }
    
    void seed( const unsigned int s ) 
    { 
        // pas de seeds aleatoires 
    }
    
    void index( const size_t id )
    {
        m_maxdim= std::max(m_maxdim, m_dim);
        
        m_dim= 0;
        m_index= id;
        
        // associe un pixel a l'indice
        int pi= id / m_spp;
        int px= pi % m_width;
        int py= pi / m_width;
        
        assert(px >= 0 && px < m_width);
        assert(py >= 0 && py < m_height);
        
        if(0)
        {
            if(px == 0 && py ==0 && (id % m_spp) == 0)
            {
                printf("[csobol] ndim %dd / %dd, tile %d, spp %d, sobol bits %d, owen bits %d\n", 
                    m_data_ndim, m_ndim, 
                    m_tile_size, m_spp,
                    m_sobol_bits, m_owen_bits);
                
                for(int i= 0; i < m_ndim; i++)
                    if(m_data_remap[i] != -1)
                        printf("sample[%d]= data[%d], ", i, m_data_remap[i]);
                    else
                        printf("sample[%d]= no data, ", i);
                printf("\n");
                
                if(0)
                {
                    char tmp[1024];
                    sprintf(tmp, "export_pixel_%05d.dat", m_spp);
                    FILE *out= fopen(tmp, "wt");
                    
                    std::vector<Float> sample(m_ndim);
                    for(int i= 0; i < m_spp; i++)
                    {
                        const uint32_t *seeds= (*m_shared_seeds)(px % m_tile_size, py % m_tile_size, 0);
                        // genere les n dimensions du sample
                        getUniformND(sample.data(), *m_shared_sobols, seeds, m_data_ndim, i, m_sobol_bits, m_owen_bits, true);    // owen sur log2(n) + 6
                        
                        for(int d= 0; d < m_data_ndim; d++)
                            fprintf(out, "%.10f ", sample[d]);
                        fprintf(out, "\n");
                    }
                    
                    fprintf(out, "#\n");
                    
                    fclose(out);
                }                
            }
        }
        
        // recupere le sample
        assert(m_data.size() == m_data_ndim);
        // permutations owen par dimension
        const uint32_t *seeds= (*m_shared_seeds)(px % m_tile_size, py % m_tile_size, 0);
        
        // genere les n dimensions du sample
        getUniformND(m_data.data(), *m_shared_sobols, seeds, m_data_ndim, id % m_spp, m_sobol_bits, m_owen_bits, true);
        
        // puis remappe les dimensions...
        m_values[0]= Float(0.5);    // centre du pixel
        m_values[1]= Float(0.5);
        m_values[2]= Float(0.5);    // centre de la lentille
        m_values[3]= Float(0.5);
        for(int i= 4; i < m_ndim; i++) 
            m_values[i]= Float(0);   // force un zero dans les dimensions non utilisees
        
        for(int i= 0; i < m_ndim; i++)
            if(m_data_remap[i] != -1)
                m_values[i]= m_data[m_data_remap[i]];
        
        // force une stratification sur les dimensions xy
        m_values[0]= (Float(px) + m_values[0]) / Float(m_width);
        m_values[1]= (Float(py) + m_values[1]) / Float(m_height);        
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);      // ne doit pas arriver, c'est index() qui fait tout le boulot...
    }

protected:
#if 0
    void getUniformND( std::vector<Float>& values,
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
        const std::vector<unsigned int>& seeds,
        const uint32_t nDims,
        const uint32_t n,
		const uint32_t nbits,
        const uint32_t owen_tree_depth = 32,
        const bool owen_permut_flag = true )
    {
        assert(values.size() == nDims);
        assert(seeds.size() == nDims);
        //~ unsigned int seed= seeds[0];
        
        uint32_t IDcode= n;		// dim 0: take n-th point as index -> into 32-bit integer IDcode
        for(int idim = 0; idim < nDims; idim++) 
        {
            IDcode= sobols[idim].getSobolInt(IDcode);	// radix-inversion + sobol
            uint32_t res_IDcode= IDcode;				// used to calculate the resulting value
            IDcode= IDcode >> (32-nbits);				// this will be used as new index for the next dimension
            if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
                res_IDcode= OwenScrambling(res_IDcode, seeds[idim], owen_tree_depth);
                //~ res_IDcode= scramble(res_IDcode, 32, seed);    //< works !!
            
            values[idim]= (Float) res_IDcode / (Float) UINT32SOBOLNORM;	// final value (double) for this dimension
        }
    }
    
    unsigned int scramble( const unsigned int x, const int dim, const unsigned int seed )
    {
        unsigned int dseed= RNG::hash(seed) ^ RNG::hash(dim);
        m_rng.seed(dseed);	        // 1 sequence par dimension, et par run
        
        // flip root, node_index == 0, implicit rng.index(0)
        unsigned int flip= m_rng.sample_range(2) << 31;
        unsigned int code= x ^ flip;        // flip MSB
        
        for(int d= 1; d < 32; d++)
        {
            unsigned int level_base= (1u << d) -1;
            unsigned int level_offset= x >> (32 - d);	// level > 0 == 2^d nodes
            unsigned int node_index= level_base + level_offset;
            
            unsigned int seed= RNG::hash(dseed) ^ RNG::hash(node_index);
            m_rng.seed(seed);          // 1 hash par noeud par dimension, et par run
            
            unsigned int flip=m_rng.sample_range(2) << (31 - d);
            code= code ^ flip;
        }
        
        return code;
    }
    
    void getUniformND_with_seed(std::vector<Float>& values,
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
        const uint32_t nDims,
        const uint32_t n,
		const uint32_t nbits,
        const uint32_t owen_tree_seed,
        const uint32_t owen_tree_depth = 32,
        const bool owen_permut_flag = true ) 
    {
        uint32_t IDcode = n;		// dim 0: take n-th point as index -> into 32-bit integer IDcode
        for (unsigned int idim = 0; idim < nDims; idim++) {
            IDcode = sobols[idim].getSobolInt(IDcode);	// radix-inversion + sobol
            uint32_t res_IDcode = IDcode;				// used to calculate the resulting value
            IDcode = IDcode >> (32-nbits);				// this will be used as new index for the next dimension
            if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
                res_IDcode = OwenScrambling(res_IDcode, owen_tree_seed, owen_tree_depth);
            values[idim] = (Float) res_IDcode / (Float) UINT32SOBOLNORM;	// final value (double) for this dimension
        }
    }	// getUniformND
#endif
    
    const BlueTile *m_shared_seeds;
    const std::vector< SobolGenerator1D<uint32_t> > *m_shared_sobols;
    
    int m_owen_bits;
    int m_sobol_bits;
    int m_tile_size;
    int m_width;
    int m_height;
    int m_spp;
};


// owen global 
struct OwenPPSampler : public Sampler
{
    OwenPPSampler( ) : Sampler(), 
        m_sobols(nullptr), m_dimensions()
    {}
    
    explicit OwenPPSampler( const std::vector<int>& d, const size_t n, const char *sobol_filename, const unsigned int s ) : Sampler(d.size(), n, s), 
        m_sobols(nullptr), m_dimensions(d)
    {
        m_sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *m_sobols);    // read sobols from file and fill appropriate structures
    }
    
    OwenPPSampler( const OwenPPSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_sobols(sampler.m_sobols), m_dimensions(sampler.m_dimensions)
    {}
    
    OwenPPSampler& operator= ( const OwenPPSampler& ) = delete;
    OwenPPSampler( const OwenPPSampler& ) = delete;
    
    ~OwenPPSampler( ) 
    {
        printf("owenpp max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("owenpp release...\n");
        delete m_sobols;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new OwenPPSampler(*this, s);
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        // recupere le sample
        assert(m_data.size());
        assert(m_data.size() == m_dimensions.size());
        // uniquement les dimensions necessaires...
        for(int d= 0; d < m_data_ndim; d++)
        {
            assert(m_dimensions[d] != -1);
            unsigned int seed= RNG48::hash(m_seed) ^ RNG48::hash(m_dimensions[d]);
            values[d]= getOwenPlus1D_with_seed(*m_sobols, m_dimensions[d], id, seed);
        }
    }
    
protected:
    std::vector< SobolGenerator1D<uint32_t> > *m_sobols;
    std::vector<int> m_dimensions;
};




struct OwenPixelSampler : public Sampler
{
    OwenPixelSampler( ) : Sampler(), 
        m_shared_sobols(nullptr), m_seeds(), m_sobol_bits(0), m_owen_bits(0)
    {}
    
    explicit OwenPixelSampler( const int ndim, const int spp, const char *sobol_filename, const int owen_bits, const unsigned int s ) : Sampler(ndim, spp, s), 
        m_shared_sobols(nullptr), m_seeds(), m_sobol_bits(0), m_owen_bits(32)
    {
        auto *sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *sobols);    // read sobols from file and fill appropriate structures
        
        m_shared_sobols= sobols;
        m_seeds.resize(m_ndim);
        
        m_sobol_bits= std::log2(spp);
        m_owen_bits= std::min(m_sobol_bits + owen_bits, 32);
    }
    
    OwenPixelSampler( const OwenPixelSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_shared_sobols(sampler.m_shared_sobols), m_seeds(sampler.m_seeds), m_sobol_bits(sampler.m_sobol_bits), m_owen_bits(sampler.m_owen_bits)
    {}
    
    OwenPixelSampler& operator= ( const OwenPixelSampler& ) = delete;
    OwenPixelSampler( const OwenPixelSampler& ) = delete;
    
    ~OwenPixelSampler( ) 
    {
        printf("owen_pixel max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("owen_pixel release...\n");
        delete m_shared_sobols;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new OwenPixelSampler(*this, s);
    }
    
    void seed( const unsigned int s ) 
    {
        m_seed= s;
        m_rng= RNG48(s);    
    }

    void pixel_seeds( const int n, const uint32_t *seeds )
    {
        m_seeds.assign(seeds, seeds+n);
    }
    
    void index( const size_t id )
    {
        m_maxdim= std::max(m_maxdim, m_dim);
        
        m_dim= 0;
        m_index= id;
        
        // recupere le sample
        assert(m_data.size());
        // uniquement les dimensions necessaires...
        for(int d= 0; d < m_data_ndim; d++)
        {
            unsigned int seed= m_seeds[d];
            m_data[d]= getOwenPlus1D_with_seed(*m_shared_sobols, d, id, seed, m_owen_bits);
        }
        
        // puis remappe les dimensions...
        m_values[0]= Float(0.5);    // centre du pixel
        m_values[1]= Float(0.5);
        m_values[2]= Float(0.5);
        m_values[3]= Float(0.5);
        for(int i= 4; i < m_ndim; i++) 
            m_values[i]= Float(0);   // force un zero dans les dimensions non utilisees
        
        for(int i= 0; i < m_ndim; i++)
            if(m_data_remap[i] != -1)
                m_values[i]= m_data[m_data_remap[i]];
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);
    }
    
protected:
    const std::vector< SobolGenerator1D<uint32_t> > *m_shared_sobols;
    std::vector<uint32_t> m_seeds;
    
    int m_sobol_bits;
    int m_owen_bits;
};


// utilise exclusivement par blue / pixel sampler, pas de stratification forcee sur le plan image. 0 < index < spp
struct CSobolPixelSampler : public Sampler
{
    CSobolPixelSampler( ) : Sampler(), 
        m_shared_sobols(nullptr), m_seeds(), m_sobol_bits(0), m_owen_bits(0)
    {}
    
    explicit CSobolPixelSampler( const int ndim, const int spp, const char *sobol_filename, const int owen_bits, const unsigned int s ) : Sampler(ndim, spp, s), 
        m_shared_sobols(nullptr), m_seeds(), m_sobol_bits(0), m_owen_bits(32)
    {
        auto *sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *sobols);    // read sobols from file and fill appropriate structures
        
        m_shared_sobols= sobols;
        m_seeds.resize(m_ndim);
        
        m_sobol_bits= std::log2(spp);
        m_owen_bits= std::min(m_sobol_bits + owen_bits, 32);
    }
    
    CSobolPixelSampler( const CSobolPixelSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_shared_sobols(sampler.m_shared_sobols), m_seeds(sampler.m_seeds), m_sobol_bits(sampler.m_sobol_bits), m_owen_bits(sampler.m_owen_bits)
    {}
    
    CSobolPixelSampler& operator= ( const CSobolPixelSampler& ) = delete;
    CSobolPixelSampler( const CSobolPixelSampler& ) = delete;
    
    ~CSobolPixelSampler( ) 
    {
        printf("csobol_pixel max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("csobol_pixel release...\n");
        delete m_shared_sobols;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new CSobolPixelSampler(*this, s);
    }
    
    void seed( const unsigned int s ) 
    {
        m_seed= s;
        m_rng= RNG48(s);    
    }

    void pixel_seeds( const int n, const uint32_t *seeds )
    {
        m_seeds.assign(seeds, seeds+n);
    }
    
    void index( const size_t id )
    {
        m_maxdim= std::max(m_maxdim, m_dim);
        
        m_dim= 0;
        m_index= id;
        
        // recupere le sample
        assert(m_data.size() == m_data_ndim);
        // genere les n dimensions du sample
        getUniformND(m_data.data(), *m_shared_sobols, m_seeds.data(), m_data_ndim, id, m_sobol_bits, m_owen_bits, true);
        
        // puis remappe les dimensions...
        m_values[0]= Float(0.5);    // centre du pixel
        m_values[1]= Float(0.5);
        m_values[2]= Float(0.5);
        m_values[3]= Float(0.5);
        for(int i= 4; i < m_ndim; i++) 
            m_values[i]= Float(0);   // force un zero dans les dimensions non utilisees
        
        for(int i= 0; i < m_ndim; i++)
            if(m_data_remap[i] != -1)
                m_values[i]= m_data[m_data_remap[i]];
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);
    }
    
protected:
    const std::vector< SobolGenerator1D<uint32_t> > *m_shared_sobols;
    std::vector<uint32_t> m_seeds;
    
    int m_sobol_bits;
    int m_owen_bits;
};


#endif
