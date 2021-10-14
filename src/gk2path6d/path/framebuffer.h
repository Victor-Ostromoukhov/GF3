
#ifndef _FRAMEBUFFER_H
#define _FRAMEBUFFER_H

#include <array>
#include <vector>
#include <algorithm>
#include <cmath>

#include "vec.h"
#include "color.h"
#include "image.h"
#include "sampler.h"


struct Framebuffer
{
    std::vector<Color> pixels;
    std::vector<Float> weights;
    std::vector<Float> zeroes;
    
    Framebuffer( const int w, const int h ) : pixels(w*h, Black()), weights(w*h, 0), zeroes(w*h, 0), m_width(w), m_height(h) {}
    
    Framebuffer() = delete;
    Framebuffer( const Framebuffer& ) = delete;
    Framebuffer& operator= ( const Framebuffer& ) = delete;
    
    int width( ) const { return m_width; }
    int height( ) const { return m_height; }
    size_t size( ) const { return pixels.size(); }
    
    Float filter( const Float u, const Float v, const int radius= 1 )
    {
        return std::max(Float(0), radius - u) * std::max(Float(0), radius - v) / (radius * radius);
    }
    
    Float distance( const Float u, const Float v )
    {
        return std::abs(u - v);
    }
    
    // stratification du plan image : chaque echantillon contribue a plusieurs pixels.
    void splat( const Float sx, const Float sy, const Color& sample, const int radius= 1 )
    {
        int offset= std::min(int(sy), m_height -1) * m_width + std::min(int(sx), m_width -1);
        assert(offset < m_width * m_height);
        
        if(std::isinf(sample.r) || std::isinf(sample.g) || std::isinf(sample.b)
        || std::isnan(sample.r) || std::isnan(sample.g) || std::isnan(sample.b))
        {
            zeroes[offset]= zeroes[offset] +1;
            return;
        }
        
        if(sample.power() == 0)
            zeroes[offset]= zeroes[offset] +1;
        
    #if 1
        // box filter
        //~ int offset= std::min(int(sy), m_height -1) * m_width + std::min(int(sx), m_width -1);
        //~ assert(offset < m_width * m_height);
        
        //~ // naive version
        //~ pixels[offset]= pixels[offset] + sample;
        //~ weights[offset]= weights[offset] + 1;
        
        // online version
        // cf https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        weights[offset]= weights[offset] + 1;
        pixels[offset]= pixels[offset] + (sample - pixels[offset]) / weights[offset];
        
    #endif
        
    #if 0
        // meme resolution, gros filtre 
        int wxmin= std::max(0, int(std::floor(sx)) - radius);
        int wxmax= std::min(m_width, int(std::floor(sx)) + radius +1);
        int wymin= std::max(0, int(std::floor(sy)) - radius);
        int wymax= std::min(m_height, int(std::floor(sy)) + radius +1);
        
        for(int y= wymin; y < wymax; y++)
        for(int x= wxmin; x < wxmax; x++)
        {
            unsigned int offset= y * m_width + x;
            Float w= filter(distance(sx, x + .5f), distance(sy, y + .5f), radius);
            if(w == 0)
                continue;
            
            pixels[offset]= pixels[offset] + w * sample;
            weights[offset]= weights[offset] + w;
        }
    #endif
    }
    
    Color operator() ( const int x, const int y ) const
    {
        unsigned int offset= y * m_width + x;
        //~ return pixels[offset] / weights[offset];
        // online
        return pixels[offset];
    }
    
    void clear( const int x, const int y ) 
    {
        unsigned int offset= y * m_width + x;
        pixels[offset]= Black();
        weights[offset]= 0;
        zeroes[offset]= 0;
    }
    
    void clear( const unsigned int offset, const Color& color= Black() ) 
    {
        //~ assert(offset < m_width * m_height);
        pixels[offset]= color;
        weights[offset]= 0;
        zeroes[offset]= 0;
    }
    
    void clear( const int x, const int y, const Color& color ) 
    {
        unsigned int offset= y * m_width + x;
        pixels[offset]= color;
        weights[offset]= 0;
        zeroes[offset]= 0;
    }
    
    void clear( const Color& color= Black() );
    
    Color operator() ( const unsigned int offset ) const
    {
        //~ assert(offset < m_width * m_height);        
        //~ return pixels[offset] / weights[offset];
        // online
        return pixels[offset];
    }
    
    Framebuffer& add( const Framebuffer& data );
    Framebuffer& add( const Framebuffer *data ) { assert(data); return add(*data); }
    
    Image flatten( ) const;
    Image flatten( const Float gamma ) const;
    Image flatten_weights( ) const;
    Image flatten_zeroes( ) const;

protected:
    int m_width;
    int m_height;
};

#endif
