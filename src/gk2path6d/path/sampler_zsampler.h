
#ifndef _Z_SAMPLER_H
#define _Z_SAMPLER_H

#include "real.h"
#include "sampler.h"

#include "../Screenspace/ZSampler/z.h"
#include "../Screenspace/ZSampler/zhash.h"
#include "../Screenspace/ZSampler/zt_sequence.h"
#include "../Screenspace/ZSampler/ztools.h"
#include "../Samplers/BlueTile.h"

struct ZHashSampler : public Sampler
{
    ZHashSampler( ) : Sampler(), 
        m_samples() 
    {}
    
    explicit ZHashSampler( const int ndim, const int width, const int height, const size_t spp, const unsigned int s ) : Sampler(ndim, width*height*spp, s), 
        m_samples(), m_width(width), m_height(height), m_spp(spp)
    {
        assert((m_ndim % 2) == 0);
        m_samples.resize(m_ndim*m_spp);
    }
    
    ZHashSampler( const ZHashSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_samples(sampler.m_samples.size()), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}
    
    ZHashSampler& operator= ( const ZHashSampler& ) = delete;
    ZHashSampler( const ZHashSampler& ) = delete;
    
    ~ZHashSampler( ) 
    {
        printf("zhash max dimension %d\n", m_maxdim);
    }
    
    void release( ) {}
    
    Sampler *clone( const unsigned int s )
    {
        return new ZHashSampler(*this, s);
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
        if((id % m_spp) == 0)
        {
            int samplesPerPixel = RoundUpPow2(m_spp);
            int resolution = RoundUpPow2(std::max(m_width, m_height));
            int log2Resolution = Log2Int(resolution);
            assert((m_data_ndim % 2) == 0);
            
            // extrait tous les samples par paire de dimensions
            for(int pair= 0; pair < m_data_ndim / 2; pair++)
                Z2DHash(px, py, 1, samplesPerPixel, log2Resolution, &m_samples[2*pair*m_spp], pair);
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
    std::vector<double> m_samples;
    int m_width;
    int m_height;
    int m_spp;
};


struct ZSampler : public Sampler
{
    ZSampler( ) : Sampler(), 
        m_ztable(nullptr), m_samples() 
    {}
    
    explicit ZSampler( const int ndim, const int width, const int height, const size_t spp, const unsigned int s ) : Sampler(ndim, width*height*spp, s), 
        m_ztable(nullptr), m_samples(), m_width(width), m_height(height), m_spp(spp)
    {
        assert((m_ndim % 2) == 0);
        m_samples.resize(m_ndim*m_spp);
        
        m_ztable= new ZTable;
        initZTable(*m_ztable, m_ndim / 2);
        // randomZTable(*m_ztable, m_ndim / 2);    
    }
    
    ZSampler( const ZSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_ztable(sampler.m_ztable), m_samples(sampler.m_samples.size()), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}
    
    ZSampler& operator= ( const ZSampler& ) = delete;
    ZSampler( const ZSampler& ) = delete;
    
    ~ZSampler( ) 
    {
        printf("z max dimension %d\n", m_maxdim);
    }
    
    void release( )
    {
        printf("z release...\n");
        delete m_ztable;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new ZSampler(*this, s);
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
        if((id % m_spp) == 0)
        {
            int samplesPerPixel = RoundUpPow2(m_spp);
            int resolution = RoundUpPow2(std::max(m_width, m_height));
            int log2Resolution = Log2Int(resolution);
            assert((m_data_ndim % 2) == 0);
            
            // extrait tous les samples par paire de dimensions
            for(int pair= 0; pair < m_data_ndim / 2; pair++)
                Z2D(*m_ztable, px, py, 1, samplesPerPixel, log2Resolution, &m_samples[2*pair*m_spp], pair);
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
        m_values[2]= Float(0.5);
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
    
    // extrait les samples
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);      // ne doit pas arriver, c'est index() qui fait tout le boulot...
    }
    
protected:
    ZTable *m_ztable;
    std::vector<double> m_samples;
    int m_width;
    int m_height;
    int m_spp;
};


struct Sampler02 : public Sampler
{
    Sampler02( ) : Sampler(), 
        m_samples(), m_shared_seeds()
    {}
    
    explicit Sampler02( const int ndim, const int width, const int height, const size_t spp, const unsigned int s ) : Sampler(ndim, width*height*spp, s), 
        m_samples(), m_shared_seeds(nullptr), m_width(width), m_height(height), m_spp(spp)
    {
        assert((ndim % 2) == 0);
        
        // pre-alloue les samples ZT
        for(int d= 0; d < ndim/2; d++)
            m_samples.push_back( std::vector<Point2f>(spp) );
        
        // seed par pixel
        int tile_size= std::min(width, height);
        auto *seeds= new BlueTile(tile_size, 1, 1);
        seeds->random_init(s);
        
        m_shared_seeds= seeds;
    }
    
    Sampler02( const Sampler02& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_samples(sampler.m_samples), m_shared_seeds(sampler.m_shared_seeds), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}
    
    Sampler02& operator= ( const Sampler02& ) = delete;
    Sampler02( const Sampler02& ) = delete;
    
    ~Sampler02( ) 
    {
        printf("sampler 02 max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        delete m_shared_seeds;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new Sampler02(*this, s);
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
        if((id % m_spp) == 0)
        {
            int tile_size= m_shared_seeds->width();
            uint32_t seed= (*m_shared_seeds)(px % tile_size, py % tile_size, 0, 0);
            ZTSequence2D(1, m_spp, m_samples, m_samples.size(), seed);
        }
        
        assert(m_data.size() == m_data_ndim);
        assert(m_samples.size() == m_data_ndim/2);
        // uniquement les dimensions necessaires...
        for(int d= 0; d < m_data_ndim/2; d++)
        {
            m_data[2*d]= m_samples[d][id % m_spp].x;
            m_data[2*d+1]= m_samples[d][id % m_spp].y;
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
    std::vector< std::vector<Point2f> > m_samples;
    const BlueTile *m_shared_seeds;
    int m_width;
    int m_height;
    int m_spp;
};

#endif

