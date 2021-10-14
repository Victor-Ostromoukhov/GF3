
#ifndef _SAMPLERND_H
#define _SAMPLERND_H

#include <cstdio>
#include <array>

#include "real.h"
#include "sampler.h"

struct SamplerNd : public Sampler
{
    SamplerNd( ) : Sampler(), 
        m_samples(nullptr), m_sets(0) 
    {}
    
    SamplerNd( const char *filename, const int ndim, const char *remap, const char *pad, const size_t n, const unsigned int seed ) : Sampler(ndim, n, seed), 
        m_samples(nullptr), m_sets(0)
    { 
        if(!read_samples(filename, remap, pad))
            exit(1);
    }
    
    SamplerNd( const SamplerNd& sampler, const unsigned int seed ) : Sampler(sampler, seed), 
        m_samples(sampler.m_samples), m_sets(sampler.m_sets)
    {}
    // tous les threads, referencent la meme copie des samples.
    // utiliser sampler::release() pour la detruire, uniquement sur main_sampler, pas sur les samplers qui ne detiennent que la reference
    
    SamplerNd& operator= ( const SamplerNd& _sampler ) = delete;
    SamplerNd( const SamplerNd& _sampler ) = delete;
    
    ~SamplerNd( ) {}
    
    Sampler *clone( const unsigned int seed )
    {
        return new SamplerNd(*this, seed);
    }
    
    void release( ) 
    {
        printf("release data samples...\n");
        delete m_samples;
    }
    
    void sample( const size_t id, std::vector<Float>& values )
    {
        assert(values.size() == m_data_ndim);
        assert(id*m_data_ndim + m_data_ndim <= m_samples->size());
        for(int i= 0; i < m_data_ndim; i++)
            values[i]= (*m_samples)[id*m_data_ndim + i];
    }
    
    int sets( ) const { return m_sets; }
    
protected:
    std::vector<Float> *m_samples;
    int m_sets;
    
    int read_sample( FILE *in, std::vector<Float>& sample  ) 
    {
        double x;
        for(int i= 0; i < int(sample.size()); i++)
        {
            if(fscanf(in, "%lf ", &x) != 1) 
                return i;
            
            sample[i]= Float(x);
        }
        
        return int(sample.size());
    }
    
    int read_sample_binary( FILE *in, std::vector<Float>& sample  ) 
    {
        double x;
        for(int i= 0; i < int(sample.size()); i++)
        {
            if(fread(&x, sizeof(double), 1, in) != 1)
                return i;
            
            sample[i]= Float(x);
        }
        
        return int(sample.size());
    }

    bool read_header_binary( FILE *in, const int ndim )
    {
        size_t header_size= 2 * sizeof(double) * ndim + sizeof(int) + sizeof(unsigned int);
        assert(header_size < 4096);
        
        unsigned char header[4096];
        return (fread(header, 1, header_size, in) == header_size);
    }
    
    bool read_samples( const char *filename, const char *remap_string, const char *pad_string )
    {
        bool binary= false;
        FILE *in= nullptr;
        
        if(std::string(filename).rfind(".bin") == std::string::npos)
        {
            in= fopen(filename, "rt");
            binary= false;
        }
        else
        {
            in= fopen(filename, "rb");
            binary= true;
        }
        
        if(!in)
        {
            printf("[error] loading samples '%s'...\n", filename);
            return false;
        }
        
        if(binary)
            printf("loading binary samples '%s'...\n", filename);
        else
            printf("loading samples '%s'...\n", filename);
        
        // remarque : le sampler lit le fichier depuis le constructeur, a ce moment la, m_ndim est egal a la dimension des samples
        // la dimension des samples est changee par l'appel a remap_dimensions(), fait plus tard dans le main
        // lors de l'utilisation de sample() il y a bien m_data_ndim dimensions (data) remappees vers m_ndim (sample)
        
        m_samples= new std::vector<Float>();
        m_samples->reserve(m_ndim*size());
        
        std::vector<Float> tmp(m_ndim, 0);
        while(true)
        {
            if(binary)
            {
                if(!read_header_binary(in, m_ndim))
                    break;
            }
            
            size_t k= 0;
            size_t last_size= m_samples->size();
            bool valid= true;
            for(; k < m_size; k++)
            {
                if(binary)
                {
                    if(read_sample_binary(in, tmp) < m_ndim) 
                        break;
                }
                else if(read_sample(in, tmp) < m_ndim)
                    break;
                
                for(int i= 0; i < m_ndim; i++)
                {
                    if(std::isnan(tmp[i]))
                        valid= false;
                    
                    if(tmp[i] < 0) tmp[i]= 0;
                    
                #ifdef FLOAT32
                    if(tmp[i] >= 1) tmp[i]= 0.99999994f;
                #else
                    if(tmp[i] >= 1) tmp[i]= 0.99999999999999989;
                #endif
                }
                
                for(int i= 0; i < m_ndim; i++) 
                    m_samples->push_back(tmp[i]);
            }
            
            if(k == 0)  // separateur ascii mal place, a la fin du fichier...
                break;
            
            if(k != m_size)
            {
                printf("[error] invalid number of points, %d points %dd != %d points, last points %d\n", int(m_samples->size() / m_ndim), m_ndim, int(m_size), int(k));
                fclose(in);
                return false;
            }
            
            if(valid)
                m_sets++;
            else
                m_samples->resize(last_size);
            
            if(!binary)
            {
                // essaye de lire un separateur
                char delim= 0;
                if(fscanf(in, "%c", &delim) != 1 || delim != '#')
                    break;
            }
        }
        fclose(in);
        
        printf("  %d sets, %d %dD points\n", m_sets, int(m_size), m_ndim);
        return true;
    }
};

#endif
