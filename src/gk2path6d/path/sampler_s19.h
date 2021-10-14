
#ifndef Sampler19_H
#define Sampler19_H

#include "sampler.h"


struct Sampler19 : public Sampler
{
    Sampler19( ) : Sampler(),
        m_dimensions(), m_optimized(nullptr), m_width(0), m_height(0), m_spp(0)
    {}
    
    explicit Sampler19( const std::vector<int>& d, const int width, const int height, const size_t spp, const unsigned int seed ); // cf sampler_s19.cpp
    
    explicit Sampler19( const Sampler19& sampler, const unsigned int seed ) : Sampler(sampler, seed),
        m_dimensions(sampler.m_dimensions), m_optimized(sampler.m_optimized), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}

    Sampler19( const Sampler19& _sampler ) = delete;
    Sampler19& operator= ( const Sampler19& _sampler ) = delete;
    
    ~Sampler19( )
    {
        printf("sampler s19 max dimension %d\n", m_maxdim);
    }
    
    void seed( const unsigned int s ) 
    { 
        // pas de seeds aleatoires 
    }
    
    Sampler *clone( const unsigned int seed )
    {
        return new Sampler19(*this, seed);
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
        
        assert(px >= 0 && px < m_width);
        assert(py >= 0 && py < m_height);
        
        // recupere le sample
        assert(m_data.size());
        assert(m_data.size() == m_dimensions.size());
        assert(m_optimized);
        // uniquement les dimensions necessaires...
        for(int d= 0; d < m_data_ndim; d++)
        {
            assert(m_dimensions[d] != -1);
            for(int d= 0; d < m_data_ndim; d++)
                m_data[d]= m_optimized(px, py, id % m_spp, m_dimensions[d]);
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
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(0);      // ne doit pas arriver, c'est index() qui fait tout le boulot...
    }

protected:
    std::vector<int> m_dimensions;    
    float (*m_optimized)(int x, int y, int i, int d );
    int m_width;
    int m_height;
    int m_spp;
};

#endif
