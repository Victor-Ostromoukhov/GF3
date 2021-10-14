
#ifndef _MORTON_SAMPLER_H
#define _MORTON_SAMPLER_H


struct MortonSampler : public Sampler
{
    MortonSampler( ) : Sampler(), 
        m_sobols(nullptr), m_width(0), m_height(0), m_spp(0) 
    {}
    
    explicit MortonSampler( const int ndim, const int width, const int height, const size_t spp, 
        const char *sobol_filename, const int sobol_shift, const unsigned int s ) : 
        Sampler(ndim, width*height*spp, s), 
        m_sobols(nullptr), m_sobol_indices(nullptr), m_morton_depth(0), m_sobol_shift(sobol_shift), m_width(width), m_height(height), m_spp(spp)
    {
        m_sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *m_sobols);    // read sobols from file and fill appropriate structures
        
        if(width != height)
        {
            printf("[error] morton sampler: not a square image %dx%d...\n", width, height);
            exit(1);
        }
        
        m_morton_depth= std::log2(std::max(width, height));
        
        // pre-calcule les indices de sobol associes aux pixels de l'image
        m_sobol_indices= new std::vector<unsigned int>();
        
    #if 1
        // morton shuffle / "blue" noise
        for(int y= 0; y < height; y++)
        for(int x= 0; x < width; x++)
            (*m_sobol_indices).push_back( shuffle_morton_tree(x, y, m_morton_depth, m_seed) );
    #else
        
        // random shuffle / white noise
        unsigned int i= 0;
        for(int y= 0; y < height; y++)
        for(int x= 0; x < width; x++)
            (*m_sobol_indices).push_back( i++ );
        
        std::default_random_engine rng(m_seed);
        std::shuffle(m_sobol_indices->begin(), m_sobol_indices->end(), rng);
    #endif
    }
    
    MortonSampler( const MortonSampler& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_sobols(sampler.m_sobols), m_sobol_indices(sampler.m_sobol_indices),
        m_morton_depth(sampler.m_morton_depth), m_sobol_shift(sampler.m_sobol_shift), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp) 
    {}
    // tous les threads, referencent la meme copie des seeds / des sobols.
    // utiliser sampler::release() pour la detruire, uniquement sur main_sampler, pas sur les samplers qui ne detiennent que la reference
    
    MortonSampler& operator= ( const MortonSampler& ) = delete;
    MortonSampler( const MortonSampler& ) = delete;
    
    ~MortonSampler( ) 
    {
        printf("morton max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("morton release...\n");
        delete m_sobols;
        delete m_sobol_indices;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new MortonSampler(*this, s);
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
        assert(m_data.size() == size_t(m_data_ndim));

        // uniquement les dimensions necessaires...
        //~ unsigned int index= shuffle_morton_tree(px, py, m_morton_depth, m_seed) * m_spp + id % m_spp;       // reclacule l'index sobol a chaque sample
        unsigned int index= (*m_sobol_indices)[py*m_width + px] * m_spp + m_sobol_shift*m_spp + id % m_spp;                           // re-utilise l'index precalcule
        // les samples d'un pixel forment une sequence d'indices de sobol, spp doit etre une puissance de 2...
        
        for(int d= 0; d < m_data_ndim; d++)
            m_data[d] = double((*m_sobols)[d].getSobolInt(index)) / double(UINT32SOBOLNORM);
        
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
    
    std::vector< SobolGenerator1D<uint32_t> > *m_sobols;
    
    std::vector<unsigned int> *m_sobol_indices;
    int m_morton_depth;
    int m_sobol_shift;
    int m_width;
    int m_height;
    int m_spp;
};


struct MortonSampler01 : public Sampler
{
    MortonSampler01( ) : Sampler(), 
        m_sobols(nullptr), m_width(0), m_height(0), m_spp(0) 
    {}
    
    explicit MortonSampler01( const int ndim, const int width, const int height, const size_t spp, 
        const char *sobol_filename, const int sobol_shift, const unsigned int s ) : 
        Sampler(ndim, width*height*spp, s), m_seeds(),
        m_sobols(nullptr), m_morton_depth(0), m_sobol_shift(sobol_shift), m_width(width), m_height(height), m_spp(spp)
    {
        if(ndim % 2)
        {
            printf("[error] morton01 sampler: %dD, odd dimensions...\n", ndim);
            exit(1);
        }
        
        m_sobols= new std::vector< SobolGenerator1D<uint32_t> >();
        loadSobolsFromFile(sobol_filename, *m_sobols);    // read sobols from file and fill appropriate structures
        
        if(width != height)
        {
            printf("[error] morton01 sampler: not a square image %dx%d...\n", width, height);
            exit(1);
        }
        
        m_morton_depth= std::log2(std::max(width, height));
        
        for(int i= 0; i < ndim / 2; i++)
            m_seeds.push_back( m_rng.sample_uint() );
    }
    
    MortonSampler01( const MortonSampler01& sampler, const unsigned int s ) : Sampler(sampler, s),
        m_sobols(sampler.m_sobols), m_seeds(sampler.m_seeds),
        m_morton_depth(sampler.m_morton_depth), m_sobol_shift(sampler.m_sobol_shift), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp) 
    {}
    // tous les threads, referencent la meme copie des seeds / des sobols.
    // utiliser sampler::release() pour la detruire, uniquement sur main_sampler, pas sur les samplers qui ne detiennent que la reference
    
    MortonSampler01& operator= ( const MortonSampler01& ) = delete;
    MortonSampler01( const MortonSampler01& ) = delete;
    
    ~MortonSampler01( ) 
    {
        printf("morton01 max dimension %d\n", m_maxdim);
    }
    
    void release( ) 
    {
        printf("morton01 release...\n");
        delete m_sobols;
    }
    
    Sampler *clone( const unsigned int s )
    {
        return new MortonSampler01(*this, s);
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
        assert(m_data.size() == size_t(m_data_ndim));
        
        // uniquement les dimensions necessaires...
        // les samples d'un pixel forment une sequence d'indices de sobol, spp doit etre une puissance de 2...
        for(int pair= 0; pair < m_data_ndim / 2; pair++)
        {
            //~ unsigned int index= shuffle_morton_tree(px, py, m_morton_depth, m_seed ^ RNG48::hash(pair)) * m_spp + m_sobol_shift*m_spp + id % m_spp;
            unsigned int index= shuffle_morton_tree(px, py, m_morton_depth, m_seeds[pair]) * (m_sobol_shift +1) * size_t(m_spp) + id % m_spp;
            //~ unsigned int index= shuffle_morton_tree(px, py, m_morton_depth, m_seeds[pair]) * m_spp + m_sobol_shift*m_spp + id % m_spp;
            // decale toute la sequence: utilise les indices [shift*main_samples .. (shift+1)*main_samples]
            
            m_data[2*pair] = double((*m_sobols)[0].getSobolInt(index)) / double(UINT32SOBOLNORM);
            m_data[2*pair+1] = double((*m_sobols)[1].getSobolInt(index)) / double(UINT32SOBOLNORM);
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
    
    std::vector< SobolGenerator1D<uint32_t> > *m_sobols;
    
    std::vector<unsigned int> m_seeds;
    int m_morton_depth;
    int m_sobol_shift;
    int m_width;
    int m_height;
    int m_spp;
};

#endif
