
#ifndef _ZPPSAMPLER_SAMPLER_H
#define _ZPPSAMPLER_SAMPLER_H

#include "real.h"
#include "sampler.h"

#include "../includes/SeedTileHelper.hpp"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"
#include "../Screenspace/ZSampler/zpp.h"
#include "../Screenspace/ZSampler/ztools.h"

struct ZPPSampler : public Sampler
{
    ZPPSampler( ) : Sampler(), 
        m_ztable(nullptr), m_sobols(nullptr), m_seeds(nullptr), m_samples(), m_width(0), m_height(0), m_spp(0)
    {}
    
    explicit ZPPSampler( const int ndim, const int width, const int height, const size_t spp, 
        const char *seed_filename, const char *sobol_filename, const int owen_bits, 
        const unsigned int s ) 
        : 
        Sampler(ndim, width*height*spp, s), 
        m_ztable(nullptr), m_sobols(nullptr), m_seeds(nullptr), m_samples(), m_width(width), m_height(height), m_spp(spp)
    {
        assert((m_ndim % 2) == 0);
        m_samples.resize(m_ndim*m_spp);
        
        m_ztable= new ZPPTable;
        initZPPTable(*m_ztable, m_ndim / 2);
        
        /*
        m_seeds= new SeedTile();
        m_seeds->loadTile(seed_filename);
        m_tile_size= m_seeds->size;
        if(m_seeds->dimension != m_ndim)
        {
            printf("[error] seed tile %dd != sampler %dd\n", m_seeds->dimension, m_ndim);
            exit(1);
        }
        */ 
        
        m_sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *m_sobols);    // read sobols from file and fill appropriate structures
        if(owen_bits == 0)
            m_owen_bits= std::log2(m_spp);
    }
    
    ZPPSampler( const ZPPSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_ztable(sampler.m_ztable), m_sobols(sampler.m_sobols), m_seeds(sampler.m_seeds), m_samples(sampler.m_samples.size()), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}
    
    ZPPSampler& operator= ( const ZPPSampler& ) = delete;
    ZPPSampler( const ZPPSampler& ) = delete;
    
    ~ZPPSampler( ) 
    {
        printf("zpp max dimension %d\n", m_maxdim);
    }
    
    void release( )
    {
        printf("zpp release...\n");
        delete m_ztable;
        delete m_sobols;
        delete m_seeds;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new ZPPSampler(*this, s);
    }
    
    // prepare tous les samples du pixel
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
        
        // ne generer les samples que lors du changement de pixel...
        // \todo a verifier, ne va pas fonctionner en multithread pour plus de 4096spp, cf decoupage de l'espace des samples par main et openmp
        if((id % m_spp) == 0)
        {
            int samplesPerPixel = RoundUpPow2(m_spp);
            int resolution = RoundUpPow2(std::max(m_width, m_height));
            int log2Resolution = Log2Int(resolution);
            assert((m_data_ndim % 2) == 0);
            
            // extrait tous les samples par paire de dimensions
            for(int pair= 0; pair < m_data_ndim / 2; pair++)
                ZPP2D(*m_ztable, px, py, 1, samplesPerPixel, log2Resolution, &m_samples[2*pair*m_spp], pair);
        }
        
        assert(m_data.size() == m_data_ndim);
        // uniquement les dimensions necessaires...
        for(int pair= 0; pair < m_data_ndim / 2; pair++)
        {
            int s= id % m_spp;
            m_data[2*pair]= m_samples[2*(pair*m_spp + s)];
            m_data[2*pair+1]= m_samples[2*(pair*m_spp + s) +1];
            // les samples sont organise par spp, puis par dimension
        }
        
        // puis remappe les dimensions...
        m_values[0]= Float(0.5);    // centre du pixel
        m_values[1]= Float(0.5);
        for(int i= 2; i < m_ndim; i++) 
            m_values[i]= Float(0);   // force un zero dans les dimensions non utilisees
            
        for(int i= 0; i < m_ndim; i++)
            if(m_data_remap[i] != -1)
                m_values[i]= m_data[m_data_remap[i]];
        
        // force une stratification sur les dimensions xy
        m_values[0]= (Float(px) + m_values[0]) / Float(m_width);
        m_values[1]= (Float(py) + m_values[1]) / Float(m_height);        
    }
    
    // extrait les samples
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);      // ne doit pas arriver, c'est index() qui fait tout le boulot...
    }
    
protected:
    ZPPTable *m_ztable;
    std::vector< SobolGenerator1D<uint32_t> > *m_sobols;
    SeedTile *m_seeds;
    
    std::vector<double> m_samples;
    
    int m_owen_bits;
    int m_tile_size;
    int m_width;
    int m_height;
    int m_spp;
};

#endif

