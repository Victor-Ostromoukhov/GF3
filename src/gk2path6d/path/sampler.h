
#ifndef  _SAMPLER_H
#define _SAMPLER_H

#include <random>

#include "real.h"
#include "image.h"
#include "image_io.h"

//~ #define WRITE_SAMPLES

// reproduit lrand48/drand48
template< typename T, T rng_a, T rng_c, T rng_m >
struct RNGT48
{
    RNGT48( ) : m_state(T(1)), m_seed(T(1)) {}
    RNGT48( const unsigned int _s ) : m_state(T(1)), m_seed(T(1)) { seed(_s); }

    RNGT48( const RNGT48& b ) = default;
    RNGT48( RNGT48&& b ) = default;
    RNGT48& operator= ( const RNGT48& b ) = default;
    RNGT48& operator= ( RNGT48&& b ) = default;

    void seed( const unsigned int _seed )
    {
        m_seed= (T(_seed) << 16) | 0x330e;      // cf srand48()
        m_state= m_seed;
    }

    void state( const T _seed )
    {
        m_seed= _seed;
        m_state= _seed;
    }

    void index( const size_t _index )
    {
        // advance to sample _index
        // cf http://www.pcg-random.org, c++ implementation
        T cur_mul= rng_a;
        T cur_add= rng_c;
        T acc_mul= T(1);
        T acc_add= T(0);
        size_t delta= _index;
        while(delta)
        {
            if(delta & 1u)
            {
                acc_mul= acc_mul * cur_mul;
                acc_add= acc_add * cur_mul + cur_add;
            }

            cur_add= cur_mul * cur_add + cur_add;
            cur_mul= cur_mul * cur_mul;
            delta= delta >> 1u;
        }

        //~ // advance current state
        //~ m_state= acc_mul * m_state + acc_add;

        // advance to sample _index
        m_state= acc_mul * m_seed + acc_add;
    }

    Float operator() ( ) { return sample(); }

    T sample_integer( )
    {
        m_state= (m_state * rng_a + rng_c) % rng_m;
        return m_state;
    }

    unsigned int sample_uint( )
    {
        T integer= sample_integer();
        return integer >> 15;   // todo verifier que T est bien 64bits... sinon return integer;
    }
    
    Float sample( )
    {
        return Float(sample_integer()) / Float(rng_m);
    }

    T sample_range( const T range )
    {
        // Efficiently Generating a Number in a Range
        // cf http://www.pcg-random.org/posts/bounded-rands.html
        T divisor= rng_m / range;
        if(divisor == 0) return 0;

        while(true)
        {
            T x= sample_integer() / divisor;
            if(x < range) return x;
        }
    }

    // c++11 random
    typedef T result_type;
    static constexpr result_type min( )  { return T(0); }
    static constexpr result_type max( )  { return T(rng_m -1); }

    // cf http://www.burtleburtle.net/bob/hash/integer.html
    static unsigned int hash( const unsigned int x )
    {
        unsigned int a= x;
        a= (a + 0x7ed55d16) + (a<<12);
        a= (a ^ 0xc761c23c) ^ (a>>19);
        a= (a + 0x165667b1) + (a<<5);
        a= (a + 0xd3a2646c) ^ (a<<9);
        a= (a + 0xfd7046c5) + (a<<3);
        a= (a ^ 0xb55a4f09) ^ (a>>16);
        
        return a;
    }
    
    static uint32_t hash3( uint32_t x )
    {
        // finalizer from murmurhash3
        x ^= x >> 16;
        x *= 0x85ebca6bu;
        x ^= x >> 13;
        x *= 0xc2b2ae35u;
        x ^= x >> 16;
        return x;
    }
    
protected:
    T m_state;
    T m_seed;
};

// cf https://en.wikipedia.org/wiki/Linear_congruential_generator
// glibc constants, 32bits
// typedef RNGT48<unsigned int, 1103515245u, 12345u, 1u << 31> RNGC;

// drand48 constants, 48bits
typedef RNGT48<uint64_t, 25214903917llu, 11llu, 1llu << 48> RNG48;

// sampler interface
struct Sampler
{
    Sampler( ) : m_ndim(0), m_size(0),
        m_dim(0), m_index(0),
        m_seed(0),
        m_maxdim(0),
        m_rng(),
        m_values(),
        m_rotation(nullptr), m_set_rotation(nullptr), m_rotation_width(0), m_rotation_height(0), m_rotation_spp(0),
        m_width(0), m_height(0), m_tile_width(0), m_tile_height(0), m_tile_spp(0),
        m_data_remap(), m_data(), m_data_ndim(0)
    {}

    explicit Sampler( const int ndim, const size_t n, const unsigned int s= 0 ) : m_ndim(ndim), m_size(n),
        m_dim(0), m_index(0),
        m_seed(s),
        m_maxdim(0),
        m_rng(s),
        m_values(ndim, 0),
        m_rotation(nullptr), m_set_rotation(nullptr), m_rotation_width(0), m_rotation_height(0), m_rotation_spp(0),
        m_width(0), m_height(0), m_tile_width(0), m_tile_height(0), m_tile_spp(0),
        m_data_remap(), m_data(), m_data_ndim(0)
    {}
    
    Sampler( const Sampler& sampler, const unsigned int _s ) : m_ndim(sampler.m_ndim), m_size(sampler.m_size),
        m_dim(0), m_index(0),
        m_seed(_s),
        m_maxdim(0),
        m_rng(_s),
        m_values(sampler.m_ndim, 0),
        m_rotation(sampler.m_rotation), m_set_rotation(sampler.m_set_rotation), m_rotation_width(sampler.m_rotation_width), m_rotation_height(sampler.m_rotation_height), m_rotation_spp(sampler.m_rotation_spp),
        m_width(sampler.m_width), m_height(sampler.m_height), m_tile_width(sampler.m_tile_width), m_tile_height(sampler.m_tile_height), m_tile_spp(sampler.m_tile_spp),
        m_data_remap(sampler.m_data_remap), m_data(sampler.m_data), m_data_ndim(sampler.m_data_ndim)
    {
        //~ printf("clone rotation(%d %d %d) %p\n", m_rotation_width, m_rotation_height, m_rotation_spp, m_rotation);
    }

    Sampler& operator= ( const Sampler& _sampler ) = delete;
    Sampler( const Sampler& _sampler ) = delete;

    virtual ~Sampler( )
    {
        //~ printf("max dimension %d\n", m_maxdim);
    }

    // les samples sont numerotes sur 3 axes :
    // dimension : par variable aleatoire
    // indice : par chemin
    // jitter / scramble / seed : par pixel

    Float sample1( ) { return sample_ndim( m_dim++, m_index ); }
    Float sample2( ) { return sample_ndim( m_dim++, m_index ); }
    Float sample3( ) { return sample_ndim( m_dim++, m_index ); }
    Float sample4( ) { return sample_ndim( m_dim++, m_index ); }
    Float sample5( ) { return sample_ndim( m_dim++, m_index ); }
    Float sample6( ) { return sample_ndim( m_dim++, m_index ); }

    // renvoie une valeur generee, ou aleatoire si la dimension est trop importante
    Float sample_ndim( const int d, const size_t id )
    {
        assert(id == m_index);

    #if 0
        if(d < m_ndim)
            return m_values[d];
        else
            return m_rng.sample();
    #else
        assert(d < m_ndim);
        return m_values[d];
    #endif
    }

    // genere toutes les dimensions
    virtual void sample( const size_t id, std::vector<Float>& values ) = 0;

    virtual Sampler *clone( const unsigned int seed ) = 0;
    virtual void release( ) {}

    size_t size( ) const { return m_size; }
    int dimensions( ) const { return m_ndim; }

    int dimension( ) const { return m_dim; }
    void dimension( const int d ) { m_dim= d; }
    size_t index( ) const { return m_index; }
    unsigned int seed( ) const { return m_seed; }

    virtual void seed( const unsigned int s )
    {
        m_seed= s;
        m_rng= RNG48(s);
        if(m_rotation)
            create_random_rotation(m_rotation_width, m_rotation_height, m_rotation_spp);
        if(m_set_rotation)
            create_set_rotation(m_rotation_width, m_rotation_height, m_rotation_spp);
    }
    
    virtual void pixel_seeds( const int n, const uint32_t *seeds ) { }  // rien !!
    // uniquement utilise par OwenPixelSampler / CSobolPixelSampler pour le blue noise brute force a posteriori, cf blue.cpp

    virtual int sets( ) const { return 1; }

    // iteration sur les samples; impl par defaut
    virtual void index( const size_t id )
    {
        m_maxdim= std::max(m_maxdim, m_dim);

        m_dim= 0;
        m_index= id;

        if(m_rotation_spp > 0)
        {
            // l'image est calculee par une seule sequence de m_size samples
            // lorsque une rotation est active, le sampler ne fournit que nspp samples
            // et la rotation associee au pixel perturbe la sequence du sampler
            
            // associe un pixel a l'indice
            int pi= id / m_rotation_spp;
            int px= pi % m_rotation_width;
            int py= pi / m_rotation_width;
            
            assert(pi >= 0 && pi < m_rotation_width * m_rotation_height);
            assert(px >= 0 && px < m_rotation_width);
            assert(py >= 0 && py < m_rotation_height);
            
            if(m_set_rotation)
                // indice du set par pixel
                m_index= (*m_set_rotation)[pi] * m_rotation_spp + id % m_rotation_spp;
            else
                // indice du sample par pixel
                m_index= id % m_rotation_spp;
            
            //~ sample(m_index, m_values);
            if(m_data.empty())
                sample(m_index, m_values);
            else
            {
                // recupere le sample
                sample(m_index, m_data);
                
                // remappe les dimensions 
                m_values[0]= Float(0.5);
                m_values[1]= Float(0.5);
                for(int i= 2; i < m_ndim; i++) 
                    m_values[i]= 0;   // force un zero dans les dimensions non utilisees
                
                for(int i= 0; i < m_ndim; i++)
                    if(m_data_remap[i] != -1)
                        m_values[i]= m_data[m_data_remap[i]];
            }
            
            if(m_rotation)
            {
                // utilise la rotation pre-generee
                assert(size_t(pi * m_ndim + m_ndim) <= m_rotation->size());
                for(int d= 0; d < m_ndim; d++)
                {
                    Float r= (*m_rotation)[pi * m_ndim + d];
                    Float u= m_values[d] + r;
                    if(u < 0) u= u +1;
                    if(u >= 1) u= u - 1;
                    assert(u >= 0 && u < 1);
                    m_values[d]= u;
                }
            }
            //~ else        // rotation == 0
            
            // force une stratification sur les dimensions xy
            m_values[0]= (Float(px) + m_values[0]) / Float(m_rotation_width);
            m_values[1]= (Float(py) + m_values[1]) / Float(m_rotation_height);
        }
        else if(m_tile_spp)
        {
            // 1 tile == tile_size samples dupliques sur l'image pour produire un total de m_size samples
            
            // associe une tuile a l'index
            size_t tile_size= size_t(m_tile_width*m_tile_height) * m_tile_spp;
            int ti= id / tile_size;
            
            m_index= id % tile_size;
            
        #if 0
            sample(m_index, m_values);
            
            // force une stratification dans la tuile sur les dimensions xy
            int tx= ti % (m_width / m_tile_width);
            int ty= ti / (m_width / m_tile_width);
            int px= tx * m_tile_width;
            int py= ty * m_tile_height;
            m_values[0]= (Float(px) + m_values[0] * m_tile_width) / Float(m_width);
            m_values[1]= (Float(py) + m_values[1] * m_tile_height) / Float(m_height);
        #endif
            
            if(m_data.empty())
                sample(m_index, m_values);
            else
            {
                // recupere le sample
                sample(m_index, m_data);
                
                // remappe les dimensions 
                for(int i= 0; i < m_ndim; i++) 
                    m_values[i]= 0;   // force un zero dans les dimensions non utilisees
                
                for(int i= 0; i < m_ndim; i++)
                    if(m_data_remap[i] != -1)
                        m_values[i]= m_data[m_data_remap[i]];
                
                // force une stratification dans la tuile sur les dimensions xy
                int tx= ti % (m_width / m_tile_width);
                int ty= ti / (m_width / m_tile_width);
                int px= tx * m_tile_width;
                int py= ty * m_tile_height;
                m_values[0]= (Float(px) + m_values[0] * m_tile_width) / Float(m_width);
                m_values[1]= (Float(py) + m_values[1] * m_tile_height) / Float(m_height);
            }
        }
        else
        {
            // genere le sample...
            //~ sample(m_index, m_values);
            
            if(m_data.empty())
                sample(m_index, m_values);
            else
            {
                // recupere le sample
                sample(m_index, m_data);
                
                // remappe les dimensions 
                m_values[0]= Float(0.5);
                m_values[1]= Float(0.5);
                for(int i= 2; i < m_ndim; i++) 
                    m_values[i]= 0;   // force un zero dans les dimensions non utilisees
                    
                for(int i= 0; i < m_ndim; i++)
                    if(m_data_remap[i] != -1)
                        m_values[i]= m_data[m_data_remap[i]];
            }
        }

    #ifdef WRITE_SAMPLES
        // stocke les valeurs
        for(int d= 0; d < m_ndim; d++)
            samples.push_back(m_values[d]);
    #endif
    }

    void remap_dimensions( const int ndim, const std::vector<int>& remap )
    {
        // dimensions des samples m_ndim / cf option -d, et option --dims
        m_data_ndim= m_ndim;
        m_data_remap= remap;
        m_data.resize(m_data_ndim);
        
        m_ndim= ndim;   // dimensions des samples apres remap...
        m_values.resize(m_ndim);
        
        printf("remapping dimensions:  data %d samples %d\n", m_data_ndim, m_ndim);
    }
    
    // force l'utilisation de la meme sequence par pixel, avec une rotation cranley-patterson aleatoire, pour decorreler les pixels.
    void create_random_rotation( const int width, const int height, const int spp, const int scale= 0 )
    {
        printf("using random rotation: %dx%d %dspp, scale %d\n", width, height, spp, scale);
        
        m_rotation_width= width;
        m_rotation_height= height;
        m_rotation_spp= spp;
        
        // recalcule le nombre de sample total
        m_size= size_t(width*height) * size_t(spp);
        
        delete m_rotation;
        m_rotation= new std::vector<Float>();
        m_rotation->reserve(m_rotation_width*m_rotation_height*m_ndim);
        
        for(int y= 0; y < m_rotation_height; y++)
        for(int x= 0; x < m_rotation_width; x++)
            for(int d= 0; d < m_ndim; d++)
            {
                if(!scale)
                    m_rotation->push_back( m_rng.sample() );
                else
                    m_rotation->push_back( 1 / Float(scale) * 2 * (m_rng.sample() - Float(0.5)) / Float(spp) );
            }
    }

    // et la rotation peut etre nulle.
    void force_null_rotation( const int width, const int height, const int spp )
    {
        printf("using null rotation: %dx%d %dspp\n", width, height, spp);
        
        m_rotation_width= width;
        m_rotation_height= height;
        m_rotation_spp= spp;
        delete m_rotation;
        m_rotation= nullptr;
        
        // recalcule le nombre de sample total
        m_size= size_t(width*height) * size_t(spp);
    }

    // et la rotation de plusieurs sets de samples
    void create_set_rotation( const int width, const int height, const int spp )
    {
        printf("using %d sets rotation...\n", sets());
        assert(sets() > 1);

        m_rotation_width= width;
        m_rotation_height= height;
        m_rotation_spp= spp;
        //~ delete m_rotation;
        //~ m_rotation= nullptr;
        delete m_set_rotation;

        //~ Image mask= read_image("data/blue.png");    // selectionne le set en fonction d'un mask blue noise
        //~ assert(mask.size());

        const int n= sets();
        m_set_rotation= new std::vector<int>();
        m_set_rotation->reserve(m_rotation_width*m_rotation_height);
        for(int y= 0; y < m_rotation_height; y++)
        for(int x= 0; x < m_rotation_width; x++)
        {
            m_set_rotation->push_back( int(m_rng.sample_range(n)) );
            //~ Float m= mask(x % mask.width(), y % mask.height()).r;
            //~ m_set_rotation->push_back(int(m * n));
        }
        
        // recalcule le nombre de sample total
        m_size= size_t(width*height) * size_t(spp);
    }

    void create_tile( const int tile_width, const int tile_height, const int width, const int height, const int spp )
    {
        printf("using %dx%d tile...\n", tile_width, tile_height);
        
        if(width % tile_width || height % tile_height)
        {
            printf("[error] incorrect tile size %dx%d, image %dx%d\n", tile_width, tile_height, width, height);
            exit(0);
        }
        
        m_width= width;
        m_height= height;
        m_tile_width= tile_width;
        m_tile_height= tile_height;
        m_tile_spp= spp;
        
        // recalcule le nombre de sample total
        printf("tile size %lld, size %lld\n", size_t(tile_width*tile_height)*spp, m_size);
        assert(size_t(tile_width*tile_height)*spp == m_size);   // tile_size samples necessaires
        m_size= size_t(width*height) * spp;
    }
    
    int write_samples( const char *filename )
    {
    #ifdef WRITE_SAMPLES
        FILE *out= fopen(filename, "wt");
        if(!out)
            return -1;
        
        for(int i= 0; i < int(m_size); i++)
        {
            for(int k= 0; k < m_ndim; k++)
                fprintf(out, "%.18lf ", m_samples[i*m_ndim + k]);
            fprintf(out, "\n");
        }
        
        fclose(out);
    #endif
        return 0;
    }

protected:
    int m_ndim;              // dimensions
    size_t m_size;              // samples
    
    int m_dim;              // next dimension to sample
    size_t m_index;             // next index to sample
    unsigned int m_seed;     // seed / scramble
    
    int m_maxdim; //
    
    RNG48 m_rng;
    
    std::vector<Float> m_values;
    
private:
    std::vector<Float> *m_rotation;
    std::vector<int> *m_set_rotation;
    int m_rotation_width, m_rotation_height;
    int m_rotation_spp;
    
    int m_width, m_height;
    int m_tile_width, m_tile_height;
    int m_tile_spp;
    
protected:
    std::vector<int> m_data_remap;
    std::vector<Float> m_data;  // m_values avant remap / permutation
    int m_data_ndim;
    
#ifdef WRITE_SAMPLES
    std::vector<Float> m_samples;
#endif
};


struct UniformSampler : public Sampler
{
    UniformSampler( ) : Sampler(),
        m_last_index(0),
        m_width(0), m_height(0), m_spp(0)
    {}
    
    explicit UniformSampler( const int ndim, const int width, const int height, const size_t spp, const unsigned int seed ) : Sampler(ndim, width*height*spp, seed),
        m_last_index(0),
        m_width(width), m_height(height), m_spp(spp)
    {}
    
    explicit UniformSampler( const UniformSampler& sampler, const unsigned int seed ) : Sampler(sampler, seed),
        m_last_index(0),
        m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}
    
    UniformSampler( const UniformSampler& _sampler ) = delete;
    UniformSampler& operator= ( const UniformSampler& _sampler ) = delete;
    
    ~UniformSampler( )
    {}
    
    Sampler *clone( const unsigned int seed )
    {
        return new UniformSampler(*this, seed);
    }
    
    // genere toutes les dimensions
    void sample( const size_t id, std::vector<Float>& values )
    {
        if(id != m_last_index +1)
            // l'image est calculee par une seule sequence aleatoire,
            // replace le generateur au bon endroit dans la sequence
            m_rng.index(id*m_ndim);
        
        // evite de le faire a chaque sample...
        m_last_index= id;
        
        for(int d= 0; d < m_ndim; d++)
            values[d]= m_rng.sample();
    }

protected:
    size_t m_last_index;
    int m_width;
    int m_height;
    int m_spp;
};


struct LHSSampler : public Sampler
{
    LHSSampler( ) : Sampler(),
        m_samples(nullptr),
        m_width(0), m_height(0), m_spp(0)
    {}
    
    explicit LHSSampler( const int ndim, const int width, const int height, const size_t spp, const unsigned int seed ) : Sampler(ndim, width*height*spp, seed),
        m_samples(nullptr), m_width(width), m_height(height), m_spp(spp)
    {
        //~ // force cranley-patterson rotation
        //~ create_random_rotation(width, height, spp);
        
        // pixel LHS
        std::vector<Float> tmp;
        tmp.reserve(m_size*m_ndim);
        
        // genere les n samples par dimension
        for(int d= 0; d < m_ndim; d++)
            for(size_t i= 0; i < m_size; i++)
                tmp.push_back( Float(i + m_rng.sample()) / Float(m_spp) );     // jitter
                //~ tmp.push_back( Float(i) / Float(m_size) );
        
        // shuffle dimension par dimension
        for(int d= 0; d < m_ndim; d++)
            std::shuffle(tmp.data() + d*m_size, tmp.data() + (d+1)*m_size, m_rng);
        
        // reconstruit les points nd / transpose tmp
        m_samples= new std::vector<Float>();
        m_samples->reserve(m_size*m_ndim);
        for(size_t i= 0; i < m_size; i++)
            for(int d= 0; d < m_ndim; d++)
                m_samples->push_back( tmp[d*m_size + i] );
    }

    explicit LHSSampler( const LHSSampler& sampler, const unsigned int seed ) : Sampler(sampler, seed),
        m_samples(sampler.m_samples), m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}

    LHSSampler( const LHSSampler& _sampler ) = delete;
    LHSSampler& operator= ( const LHSSampler& _sampler ) = delete;

    ~LHSSampler( )
    {
        printf("LHS max dimension %d\n", m_maxdim);
    }
    Sampler *clone( const unsigned int seed )
    {
        return new LHSSampler(*this, seed);
    }

    void release( )
    {
        printf("release LHS samples...\n");
        delete m_samples;
    }

    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(id*m_ndim + m_ndim <= m_samples->size());
        for(int i= 0; i < m_ndim; i++)
            values[i]= (*m_samples)[id*m_ndim + i];
    }

protected:
    std::vector<Float> *m_samples;
    int m_width;
    int m_height;
    int m_spp;
};


struct RNGSampler : public Sampler
{
    RNGSampler( ) : Sampler(),
        m_last_index(0),
        m_width(0), m_height(0), m_spp(0)
    {}

    explicit RNGSampler( const int ndim, const int width, const int height, const size_t spp, const unsigned int seed ) : Sampler(ndim, width*height*spp, seed),
        m_last_index(0),
        m_width(width), m_height(height), m_spp(spp)
    {}

    explicit RNGSampler( const RNGSampler& sampler, const unsigned int seed ) : Sampler(sampler, seed),
        m_last_index(0),
        m_width(sampler.m_width), m_height(sampler.m_height), m_spp(sampler.m_spp)
    {}

    RNGSampler( const RNGSampler& _sampler ) = delete;
    RNGSampler& operator= ( const RNGSampler& _sampler ) = delete;

    ~RNGSampler( )
    {
        printf("rng max dimension %d, %dx%d %dspp\n", m_maxdim, m_width, m_height, m_spp);
    }
    
    Sampler *clone( const unsigned int seed )
    {
        return new RNGSampler(*this, seed);
    }

    // genere toutes les dimensions
    void sample( const size_t id, std::vector<Float>& values )
    {
        if(id != m_last_index +1)
            // l'image est calculee par une seule sequence aleatoire,
            // replace le generateur au bon endroit dans la sequence
            m_rng.index(id*m_ndim);

        // evite de le faire a chaque sample...
        m_last_index= id;

        for(int d= 0; d < m_ndim; d++)
            values[d]= m_rng.sample();

        // force une stratification sur les dimensions xy
        // associe un pixel a l'indice
        int pi= id / m_spp;
        int px= pi % m_width;
        int py= pi / m_width;
        m_values[0]= (Float(px) + m_values[0]) / Float(m_width);
        m_values[1]= (Float(py) + m_values[1]) / Float(m_height);
        // \todo ne fonctionne pas avec les options remap_dims... les dimensions 0 et 1 des samples ne sont plus forcement utilisees pour sampler le plan image...
    }

protected:
    size_t m_last_index;
    int m_width;
    int m_height;
    int m_spp;
};


struct DebugSampler : public Sampler
{
    DebugSampler( ) : Sampler(), m_last_index(0) {}

    explicit DebugSampler( const int ndim, const int width, const int height, const int spp, const unsigned int s ) : Sampler(ndim, width*height*spp, s),
        m_last_index(0),
        m_grid(0), m_grid_size(0), m_spp(spp)
    {
        m_grid= std::sqrt(width);
        m_grid_size= m_grid*m_grid;
        if(m_grid_size != width || width != height)
        {
            printf("[error] debug %dd strata... wrong dimensions %dx%d\n", m_ndim, width, height);
            exit(1);
        }

        if(m_ndim == 4)
            printf("debug 4d %dx%d strata\n", m_grid, m_grid);
        else if(m_ndim == 6)
            printf("debug 6d %dx%d strata\n", m_grid, m_grid);
        else
        {
            printf("[error] debug %dd strata... not 4d or 6d\n", m_ndim);
            exit(1);
        }
    }

    DebugSampler( const DebugSampler& sampler, const unsigned s ) : Sampler(sampler, s),
        m_last_index(0),
        m_grid(sampler.m_grid), m_grid_size(sampler.m_grid_size), m_spp(sampler.m_spp)
    {}

    ~DebugSampler( ) { printf("debug max dimension %d\n", m_maxdim); }

    Sampler *clone( const unsigned int s ) { return new DebugSampler(*this, s); }

    void sample( const size_t index, std::vector<Float>& values )
    {
        //~ if(index != m_last_index +1)
            //~ m_rng.index(index);
        //~ m_last_index= index;

        size_t id= index / m_spp;
        assert(id < size_t(m_grid_size*m_grid_size));

        int grid_id= id / m_grid_size;
        int grid_px= grid_id % m_grid;
        int grid_py= grid_id / m_grid;
        assert(grid_py < m_grid);

        int cell_id= id % m_grid_size;
        int cell_px= cell_id % m_grid;
        int cell_py= cell_id / m_grid;
        assert(cell_py < m_grid);

        int px= grid_px*m_grid + cell_px;
        int py= grid_py*m_grid + cell_py;

        //~ printf("id %lu offset %d, grid %d, %d, cell %d, %d\n", id, py * m_grid_size + px, grid_px, grid_py, cell_px, cell_py);

        //~ values[0]= Float(px + 0.5) / Float(m_grid_size);
        //~ values[1]= Float(py + 0.5) / Float(m_grid_size);
        values[0]= Float(px + 0.5);
        values[1]= Float(py + 0.5);
        assert(values[0] < m_grid_size);
        assert(values[1] < m_grid_size);

        //~ values[2]= Float(grid_px + m_rng.sample()) / Float(m_grid);
        //~ values[3]= Float(grid_py + m_rng.sample()) / Float(m_grid);
        //~ values[2]= Float(grid_px + 0.5) / Float(m_grid);
        //~ values[3]= Float(grid_py + 0.5) / Float(m_grid);
        values[2]= Float(grid_px) / Float(m_grid);
        values[3]= Float(grid_py) / Float(m_grid);

        if(int(values.size()) >= 6)
        {
            //~ values[4]= Float(cell_px + m_rng.sample()) / Float(m_grid);
            //~ values[5]= Float(cell_py + m_rng.sample()) / Float(m_grid);
            //~ values[4]= Float(cell_px + 0.5) / Float(m_grid);
            //~ values[5]= Float(cell_py + 0.5) / Float(m_grid);
            values[4]= Float(cell_px) / Float(m_grid);
            values[5]= Float(cell_py) / Float(m_grid);
        }

        //~ if(cell_id == 0) printf("grid %d %d %f %f\n", grid_px, grid_py, values[2], values[3]);
        
        // \todo ne fonctionne pas avec les options remap_dims...
    }

protected:
    size_t m_last_index;

public:
    int m_grid;
    int m_grid_size;
    int m_spp;
};

#endif
