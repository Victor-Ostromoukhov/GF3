
#ifndef _PATH_TEXTURE_H
#define _PATH_TEXTURE_H

#include <vector>

#include "mesh.h"
#include "image_io.h"


struct Texture
{
    std::vector<ImageData> mipmaps;
    int base_level;
    
    Texture( ) : base_level(0) {}
    Texture( const ImageData& data );
    void build_mipmaps( const int base );
    
    int width( const Float lod= 0 ) const;
    int height( const Float lod= 0 ) const;
    int size( const Float lod= 0 ) const;
    Color sample( const Float u, const Float v, const Float lod ) const;
};


struct TextureLib
{
    std::vector<Texture> textures;
    Texture default_texture;
    bool opaque;
    
    TextureLib( );
    int insert( const ImageData& data );
    int size( ) const { return int(textures.size()); }
};


//! charge les textures associees a un ensemble de matieres, sans depasser une limite de taille, 3Go par defaut.
TextureLib read_textures( std::vector<Material>& materials, const size_t max_size= 3u*1024u*1024u*1024u );

//! detruit les textures.
void release_textures( TextureLib& lib );

#endif
