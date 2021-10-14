
#ifndef _IMAGE_H
#define _IMAGE_H

#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>

#include "color.h"


//! \addtogroup image utilitaires pour manipuler des images
///@{

//! \file
//! manipulation simplifiee d'images

//! representation d'une image.
class Image
{
protected:
    std::vector<Color32> m_data;
    int m_width;
    int m_height;

    unsigned int offset (const int x, const int y ) const
    {
        return x + y * m_width;
    }
    
public:
    Image( ) : m_data(), m_width(0), m_height(0) {}
    Image( const int w, const int h, const Color32& color= Black() ) : m_data(w*h, color), m_width(w), m_height(h) {}
    
    Image( const Image& ) = default;
    Image& operator= ( const Image& ) = default;
    
    /*! renvoie une reference sur la couleur d'un pixel de l'image.
    permet de modifier et/ou de connaitre la couleur d'un pixel :
    \code
    Image image(512, 512);
    
    image(10, 10)= make_red();      // le pixel (10, 10) devient rouge
    image(0, 0)= image(10, 10);     // le pixel (0, 0) recupere la couleur du pixel (10, 10)
    \endcode
    */
    Color32& operator() ( const int x, const int y )
    {
        return m_data[offset(x, y)];
    }
    
    //! renvoie la couleur d'un pixel de l'image (image non modifiable).
    Color32 operator() ( const int x, const int y ) const
    {
        return m_data[offset(x, y)];
    }
    
    Color32& operator() ( const unsigned int offset )
    {
        return m_data[offset];
    }
    
    //! renvoie la couleur d'un pixel de l'image (image non modifiable).
    Color32 operator() ( const unsigned int offset ) const
    {
        return m_data[offset];
    }
    
    //! renvoie la couleur interpolee a la position (x, y).
    Color32 sample( const Float x, const Float y ) const
    {
        // interpolation bilineaire 
        Float u= x - std::floor(x);
        Float v= y - std::floor(y);
        int ix= std::floor(x);
        int iy= std::floor(y);
        int ix1= std::min(ix+1, m_width -1);
        int iy1= std::min(iy+1, m_height -1);
        
        return (*this)(ix, iy)  * float((1 - u) * (1 - v))
            + (*this)(ix1, iy)  * float(u       * (1 - v))
            + (*this)(ix, iy1)  * float((1 - u) * v)
            + (*this)(ix1, iy1) * float(u       * v);
    }
    
    //! renvoie un pointeur sur le stockage des couleurs des pixels.
    const void * buffer( ) const
    {
        assert(!m_data.empty());
        return &m_data.front();
    }
    
    //! renvoie la largeur de l'image.
    int width( ) const { return m_width; }
    //! renvoie la hauteur de l'image.
    int height( ) const { return m_height; }
    //! renvoie le nombre de pixels de l'image.
    std::size_t size( ) const { return m_width * m_height; }
    
    /*! sentinelle pour la gestion d'erreur lors du chargement d'un fichier.
    exemple :
    \code
    Image image= read_image("debug.png");
    if(image == Image::error())
        return "erreur de chargement";
    \endcode
    */
    static Image& error( )
    {
        static Image image;
        return image;
    }
    
    //! comparaison avec la sentinelle. \code if(image == Image::error()) { ... } \endcode
    bool operator== ( const Image& im ) const
    {
        // renvoie vrai si im ou l'objet est la sentinelle
        return (this == &im);
    }
};

///@}
#endif
