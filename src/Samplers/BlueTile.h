
#ifndef BLUE_TILE_H
#define BLUE_TILE_H

#include <cassert>
#include <vector>
#include <random>
#include <fstream>

#include "Random.h"

struct BlueTile
{
    BlueTile() = default;
    
    BlueTile( const int size, const int ndim, const int nspp ) : m_tile(size*size*ndim*nspp), m_size(size), m_ndim(ndim), m_nspp(nspp) {}
    
    uint32_t& operator() ( const size_t id ) { return m_tile[id]; }
    uint32_t operator() ( const size_t id ) const { return m_tile[id]; }
    
    uint32_t& operator() ( const int x, const int y, const int d, const int s ) { return m_tile[offset(x, y, d, s)]; }
    uint32_t operator() ( const int x, const int y, const int d, const int s ) const  { return m_tile[offset(x, y, d, s)]; }
    
    uint32_t *operator() ( const int x, const int y, const int s ) { return m_tile.data() + offset(x, y, 0, s); }   // renvoie le pointeur sur les seeds du pixel(x, y) pour n spp
    const uint32_t *operator() ( const int x, const int y, const int s ) const { return m_tile.data() + offset(x, y, 0, s); }
    
    int width( ) const { return m_size; }
    int height( ) const { return m_size; }
    int dimensions( ) const { return m_ndim; }
    int samples() const { return m_nspp; }
    
    size_t size( ) const { return m_tile.size(); }  // renvoie nombre total de valeurs dans la tuile
    
    void random_init( const unsigned int seed= 1 )
    {
        const size_t n= m_size*m_size*m_ndim*m_nspp;
        m_tile.resize(n);
        
        RNG rng(seed);
        for(size_t i= 0; i < n; i++)
            (*this)(i)= rng.sample();
    }
    
    template< typename RNG >
    void random_init( RNG& rng )    // std::random c++ rng
    {
        const size_t n= m_size*m_size*m_ndim*m_nspp;
        m_tile.resize(n);
        
        std::uniform_int_distribution<uint32_t> u(0);   // (0, max uint32_t implicite) UINT_MAX
        for(size_t i= 0; i < n; i++)
            (*this)(i)= u(rng);
    }
    
    void read( const std::string& filename )
    {
        m_tile.clear();
        std::ifstream ifs(filename);
        if(!ifs.is_open())
        {
            std::cerr << "[error] loading tilemap '" << filename << "'...\n";
            exit(1);
        }

        std::string magic;
        ifs >> magic;
        if(magic == "#bluetile")
        {
            ifs >> m_size;
            ifs >> m_ndim;
            ifs >> m_nspp;
            
            const size_t n= m_size*m_size*m_ndim*m_nspp;
            m_tile.resize(n);
            
            std::cerr << "loading tilemap '" << filename << "' (size= "<< m_size <<", dimension= "<<  m_ndim << ", samples= " << m_nspp << ")\n";
            
            for(size_t i=0; i < n; i++)
            {
                uint32_t v;
                ifs >> v;
                (*this)(i)= v;
            }
        }
        else
        {
            std::cerr << "[error] loading tilemap '"<< filename << "' wrong format...\n";
            exit(1);
        }
        
        ifs.close();
    }
    
    void write( const std::string& filename ) const
    {
        std::ofstream ofs(filename,  std::ios::out);
        ofs << "#bluetile" << std::endl;
        ofs << m_size << " " << m_ndim << " " << m_nspp << std::endl;
        
        size_t i= 0;
        for(int y= 0; y < m_size; y++)
        for(int x= 0; x < m_size; x++)
        {
            // 1 ligne par pixel
            for(int s= 0; s < m_nspp; s++)
            for(int d= 0; d < m_ndim; d++, i++)
            {
                assert(offset(x, y, d, s) == i);    // verifie l'indexation
                ofs << (*this)(i) << " ";
            }
            ofs << std::endl;
        }
        
        ofs.close();
    }
    
    size_t offset( const int x, const int y, const int d, const int n ) const
    {
        assert(x < m_size);
        assert(y < m_size);
        assert(d < m_ndim);
        assert(n < m_nspp);
        return size_t(y*m_size + x) * (m_ndim*m_nspp) + n*m_ndim + d;
        // les seeds des dimensions d'un pixel sont en sequence...
    }

protected:
    std::vector<uint32_t> m_tile;
    int m_size;
    int m_ndim;
    int m_nspp;
};

#endif
