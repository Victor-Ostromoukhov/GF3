
#include "framebuffer.h"


Framebuffer& Framebuffer::add( const Framebuffer& data )
{
    unsigned int offset= 0;
    for(int y= 0; y < m_height; y++)
    for(int x= 0; x < m_width; x++, offset++)
    {
        Float w= weights[offset] + data.weights[offset];
        if(w > 0)
        {
            pixels[offset]= weights[offset] /  w * pixels[offset] + data.weights[offset] / w * data.pixels[offset];
            weights[offset]= w;
        }
        zeroes[offset]= zeroes[offset] + data.zeroes[offset];
    }
    
    return *this;
}

void Framebuffer::clear( const Color& color )
{
    unsigned int offset= 0;
    for(int y= 0; y < m_height; y++)
    for(int x= 0; x < m_width; x++, offset++)
    {
        pixels[offset]= color;
        weights[offset]= 0;
    }
}

Image Framebuffer::flatten( ) const
{
    Image im(m_width, m_height);
    
    unsigned int offset= 0;
    for(int y= 0; y < m_height; y++)
    for(int x= 0; x < m_width; x++, offset++)
        //~ // naive version
        //~ im(x, y)= Color32(pixels[offset] / weights[offset], 1);
        // online version
        im(x, y)= Color32(pixels[offset], 1);

    return im;
}

Image Framebuffer::flatten_weights( ) const
{
    Image im(m_width, m_height);
    
    unsigned int offset= 0;
    for(int y= 0; y < m_height; y++)
    for(int x= 0; x < m_width; x++, offset++)
        if(weights[offset] > 0)
            //~ im(x, y)= Color32(weights[offset] / 255);
            im(x, y)= Color32(weights[offset]);
        else
            im(x, y)= Red();

    return im;
}

Image Framebuffer::flatten_zeroes( ) const
{
    Image im(m_width, m_height);
    
    unsigned int offset= 0;
    for(int y= 0; y < m_height; y++)
    for(int x= 0; x < m_width; x++, offset++)
        if(weights[offset] > 0)
            im(x, y)= Color32(zeroes[offset] / weights[offset]);
            //~ im(x, y)= Color32(zeroes[offset]);
    
    return im;
}

template< typename T >
ColorT<T> pow( const ColorT<T>& x, const T y )
{
    return ColorT<T>(std::pow(x.r, y), std::pow(x.g, y), std::pow(x.b, y), 1);
}

Image Framebuffer::flatten( const Float gamma ) const
{
    Image im(m_width, m_height);
    
    unsigned int offset= 0;
    for(int y= 0; y < m_height; y++)
    for(int x= 0; x < m_width; x++, offset++)
        //~ // naive version
        //~ im(x, y)= Color32(pow(pixels[offset] / weights[offset], 1 / gamma), 1);
        // online version
        im(x, y)= Color32(pow(pixels[offset], Float(1) / gamma), 1);

    return im;
}
