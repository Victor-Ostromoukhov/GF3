
#ifndef _SOBOL_SAMPLER_H
#define _SOBOL_SAMPLER_H 

#include "real.h"
#include "sampler.h"
#include "sobol_spack.h"


struct SobolSampler : public Sampler
{
    SobolSampler( ) : Sampler() {}
    explicit SobolSampler( const int ndim, const size_t n, const unsigned int scramble ) : Sampler(ndim, n, scramble) {}
    
    SobolSampler( const SobolSampler& sampler, const unsigned int scramble ) : Sampler(sampler, scramble) {}
    
    ~SobolSampler( ) 
    { 
        //~ printf("sobol max dimension %d, scramble %08x\n", m_maxdim, m_seed); 
    }
    
    SobolSampler& operator= ( const SobolSampler& ) = delete;
    SobolSampler( const SobolSampler& ) = delete;
    
    Sampler *clone( const unsigned int seed )
    {
        return new SobolSampler(*this, seed);
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        for(int d= 0; d < m_ndim; d++)
            //~ values[d]= sobol_spack::sample(id, d, m_seed);
            values[d]= sobol_spack::sample(id, d, 0);   // deterministic sequence...
    }
};


struct OwenSampler : public Sampler
{
    OwenSampler( ) : Sampler(), 
        m_seeds(nullptr), m_scramble_rng(0), m_width(0), m_height(0), m_spp(0) 
    {}
    
    explicit OwenSampler( const int ndim, const int width, const int height, const size_t spp, const unsigned int s ) : Sampler(ndim, width*height*spp, s), 
        m_seeds(nullptr), m_scramble_rng(0), m_width(width), m_height(height), m_spp(spp)
    {
        seed(s);
    }
    
    OwenSampler( const OwenSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_seeds(sampler.m_seeds), m_scramble_rng(0), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp) 
    {}
    // tous les threads, referencent la meme copie des seeds.
    // utiliser sampler::release() pour la detruire, uniquement sur main_sampler, pas sur les samplers qui ne detiennent que la reference
    
    ~OwenSampler( ) 
    {
        printf("owen max dimension %d, scramble %08x\n", m_maxdim, m_seed);
    }
    
    OwenSampler& operator= ( const OwenSampler& ) = delete;
    OwenSampler( const OwenSampler& ) = delete;
    
    void release( ) 
    {
        printf("owen release seeds...\n");
        delete m_seeds;
    }
    
    void seed( const unsigned int s )
    {
        m_seed= s;
        m_rng= RNG48(s);
        
        delete m_seeds;
        m_seeds= new std::vector<unsigned int>(m_width*m_height);
        // seeds deterministes... 
        m_scramble_rng.seed(s);
        for(int y= 0; y < m_height; y++)
        for(int x= 0; x < m_width; x++)
            // 1 seed par pixel
            (*m_seeds)[x + y*m_width]= RNG48::hash(m_scramble_rng.sample_uint());
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new OwenSampler(*this, s);
    }

    virtual void index( const size_t id )
    {
        m_maxdim= std::max(m_maxdim, m_dim);
        
        m_dim= 0;
        m_index= id;
        
        // associe un pixel a l'indice
        int pi= id / m_spp;
        int px= pi % m_width;
        int py= pi / m_width;
        
        assert(pi >= 0 && pi < m_width * m_height);
        assert(px >= 0 && px < m_width);
        assert(py >= 0 && py < m_height);
        
        // recupere le seed du pixel
        unsigned int seed= (*m_seeds)[px + py * m_width];
        
        // recupere le sample
        if(m_data.empty())
        {
            // toutes les dimensions
            for(int d= 0; d < m_ndim; d++)
            {
                unsigned int x= sobol_spack::sample_binary(id % m_spp, d, 0);       // sequence par pixel + permutations par pixel
                unsigned int y= scramble_fast(x, d, seed);
                
                m_values[d]=  Float(y) / Float(1llu << 32);
            }
        }
        else
        {
            // uniquement les dimensions necessaires...
            for(int d= 0; d < m_data_ndim; d++)
            {
                unsigned int x= sobol_spack::sample_binary(id % m_spp, d, 0);       // sequence par pixel + permutations par pixel
                unsigned int y= scramble_fast(x, d, seed);
                
                m_data[d]=  Float(y) / Float(1llu << 32);
            }
            
            // puis remappe les dimensions...
            m_values[0]= Float(0.5);    // centre du pixel
            m_values[1]= Float(0.5);
            for(int i= 2; i < m_ndim; i++) 
                m_values[i]= Float(0);   // force un zero dans les dimensions non utilisees
                
            for(int i= 0; i < m_ndim; i++)
                if(m_data_remap[i] != -1)
                    m_values[i]= m_data[m_data_remap[i]];
        }
        
        // force une stratification sur les dimensions xy
        m_values[0]= (Float(px) + m_values[0]) / Float(m_width);
        m_values[1]= (Float(py) + m_values[1]) / Float(m_height);        
    }  

    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);      // ne doit pas arriver, c'est index() qui fait tout le boulot...
    }
    
protected:
    unsigned int scramble64( const unsigned int x, const int dim, const unsigned int seed )
    {
        //~ unsigned int dseed= RNG48::hash(m_seed) ^ RNG48::hash(dim);
        unsigned int dseed= RNG48::hash(seed) ^ RNG48::hash(dim);
        m_scramble_rng.seed(dseed);	        // 1 sequence par dimension, et par run
        
        unsigned int code= x;
        //~ unsigned int code= 0;
        for(int d= 0; d < 32; d++)
        {
            unsigned int level_base= (1u << d) -1;
            unsigned int level_offset= int64_t(x) >> (32 - d);	// level 0 == 1 node
            unsigned int node_index= level_base + level_offset;
            
            m_scramble_rng.index(node_index);	// 1 nombre aleatoire par noeud
            
            unsigned int flip=m_scramble_rng.sample_range(2) << (31 - d);
            //~ unsigned int bit= x & (1u << (31 - d));
            //~ code= code | (bit ^ flip);
            code= code ^ flip;
        }
        
        return code;
    }
    
    unsigned int scramble_fast64( const unsigned int x, const int dim, const unsigned int seed )
    {
        //~ unsigned int dseed= RNG48::hash(m_seed) ^ RNG48::hash(dim);
        unsigned int dseed= RNG48::hash(seed) ^ RNG48::hash(dim);
        
        unsigned int code= x;
        //~ unsigned int code= 0;
        for(int d= 0; d < 32; d++)
        {
            unsigned int level_base= (1u << d) -1;
            unsigned int level_offset= int64_t(x) >> (32 - d);	// level 0 == 1 node
            unsigned int node_index= level_base + level_offset;
            
            unsigned int seed= dseed ^ RNG48::hash(node_index);
            m_scramble_rng.seed(seed);          // 1 hash par noeud par dimension, et par run
            
            unsigned int flip= m_scramble_rng.sample_range(2) << (31 - d);
            //~ unsigned int bit= x & (1u << (31 - d));
            //~ code= code | (bit ^ flip);
            code= code ^ flip;
        }
        
        return code;
    }
    
    unsigned int scramble( const unsigned int x, const int dim, const unsigned int seed )
    {
        //~ unsigned int dseed= RNG48::hash(m_seed) ^ RNG48::hash(dim);
        unsigned int dseed= RNG48::hash(seed) ^ RNG48::hash(dim);
        m_scramble_rng.seed(dseed);	        // 1 sequence par dimension, et par run
        
        // flip root, node_index == 0, implicit rng.index(0)
        unsigned int flip= m_scramble_rng.sample_range(2) << 31;
        unsigned int code= x ^ flip;        // flip MSB
        
        for(int d= 1; d < 32; d++)
        {
            unsigned int level_base= (1u << d) -1;
            unsigned int level_offset= x >> (32 - d);	// level > 0 == 2^d nodes
            unsigned int node_index= level_base + level_offset;
            
            m_scramble_rng.index(node_index);	// 1 nombre aleatoire par noeud
            
            unsigned int flip=m_scramble_rng.sample_range(2) << (31 - d);
            code= code ^ flip;
        }
        
        return code;
    }
    
    unsigned int scramble_fast( const unsigned int x, const int dim, const unsigned int seed )
    {
        //~ unsigned int dseed= RNG48::hash(m_seed) ^ RNG48::hash(dim);
        unsigned int dseed= RNG48::hash(seed) ^ RNG48::hash(dim);
        m_scramble_rng.seed(dseed);	        // 1 sequence par dimension, et par run
        
        // flip root, node_index == 0, implicit rng.index(0)
        unsigned int flip= m_scramble_rng.sample_range(2) << 31;
        unsigned int code= x ^ flip;        // flip MSB
        
        for(int d= 1; d < 32; d++)
        {
            unsigned int level_base= (1u << d) -1;
            unsigned int level_offset= x >> (32 - d);	// level > 0 == 2^d nodes
            unsigned int node_index= level_base + level_offset;
            
            unsigned int seed= dseed ^ RNG48::hash(node_index);
            m_scramble_rng.seed(seed);          // 1 hash par noeud par dimension, et par run
            
            unsigned int flip=m_scramble_rng.sample_range(2) << (31 - d);
            code= code ^ flip;
        }
        
        return code;
    }

    
protected:
    std::vector<unsigned int> *m_seeds;
    RNG48 m_scramble_rng;
    int m_width;
    int m_height;
    int m_spp;
};

#endif
