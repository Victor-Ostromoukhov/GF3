
#ifndef _IMAGE_IO_H
#define _IMAGE_IO_H

#include <string>

#include "image.h"


//! \addtogroup image utilitaires pour manipuler des images
///@{

//! \file
//! manipulation directe d'images

//! charge une image a partir d'un fichier. renvoie Image::error() en cas d'echec. a detruire avec image::release( ).
//! \param filemane nom de l'image a charger
Image read_image( const char *filename );

//! enregistre une image dans un fichier png.
int write_image( const Image& image, const char *filename );


//! stockage temporaire des donnees d'une image.
struct ImageData
{
    ImageData( ) : data(), width(0), height(0), channels(0), size(0) {}
    ImageData( const int w, const int h, const int c, const int s= 1 ) : data(w*h*c*s), width(w), height(h), channels(c), size(s) {}
    
    void release( ) { std::vector<unsigned char>().swap(data); }        // detruit le vecteur / remplace par un vecteur vide.
    
    std::size_t offset( const int x, const int y ) const { return y * width * channels * size + x * channels * size; }
    const void *buffer( ) const { return &data.front(); }
    void *buffer( ) { return &data.front(); }
    
    std::vector<unsigned char> data;
    
    int width;
    int height;
    int channels;
    int size;
};

//! charge les donnees d'un fichier png. renvoie une image initialisee par defaut en cas d'echec.
ImageData read_image_data( const char *filename );

//! enregistre des donnees dans un fichier png.
int write_image_data( ImageData& image, const char *filename );

std::string path_name( const std::string& filename );
std::string base_name( const std::string& filename );

//! remplace l'extension d'un fichier
std::string change_extension( const std::string& filename, const std::string& ext );
//! ajoute un prefixe au fichier
std::string change_prefix( const std::string& prefix, const std::string& filename );

///@}

#endif
