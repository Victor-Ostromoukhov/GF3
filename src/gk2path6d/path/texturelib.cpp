
#include <cassert>

#include "texturelib.h"


static
int miplevel_size( const ImageData& image, const int lod )
{
    int w= std::max(1, image.width / (1<<lod));
    int h= std::max(1, image.height / (1<<lod));
    return image.channels * image.size * w * h;
}

static inline
int offset( const int width, const int stride, const int x, const int y, const int c )
{
    return (y * width + x) * stride + c;    // * size...
}

static
ImageData mipmap_resize( const ImageData& image )
{
    assert(image.size == 1);
    int w= std::max(1, image.width / 2);
    int h= std::max(1, image.height / 2);
    int stride= image.channels * image.size;
    
    ImageData level= ImageData(w, h, image.channels, image.size);
    
    for(int y= 0; y < h; y++)
    for(int x= 0; x < w; x++)
    for(int i= 0; i < image.channels; i++)
    {
        int m= 0;
        m= m + image.data[offset(image.width, stride, 2*x, 2*y, i)];
        m= m + image.data[offset(image.width, stride, 2*x +1, 2*y, i)];
        m= m + image.data[offset(image.width, stride, 2*x, 2*y +1, i)];
        m= m + image.data[offset(image.width, stride, 2*x +1, 2*y +1, i)];
        
        level.data[offset(w, stride, x, y, i)]= m / 4;
    }
    
    return level;
}


Texture::Texture( const ImageData& data ) : base_level(0)
{
    mipmaps.push_back(data);
}

void Texture::build_mipmaps( const int base )
{
    assert(mipmaps[0].data.size() > 0);
    int w= mipmaps[0].width;
    int h= mipmaps[0].height;
    
    // construit les mipmaps 
    for(int level= 1; w > 1 || h > 1; level++)
    {
        mipmaps.push_back(mipmap_resize(mipmaps[level -1]));
        w= mipmaps.back().width;
        h= mipmaps.back().height;
    }
    
    // ne conserve que les niveaux >= base
    base_level= base;
    for(int i= 0; i < base; i++)
    {
        mipmaps[i].release();
        mipmaps[i].width= 0;
        mipmaps[i].height= 0;
    }
}

int Texture::width( const Float lod ) const
{
    assert(mipmaps.size() > 0);
    
    int level= lod;
    if(level < base_level) level= base_level;
    if(level >= int(mipmaps.size())) level= int(mipmaps.size()) -1;
    
    assert(level < int(mipmaps.size()));
    const ImageData& image= mipmaps[level];
    return image.width;
}

int Texture::height( const Float lod ) const
{
    assert(mipmaps.size() > 0);
    
    int level= lod;
    if(level < base_level) level= base_level;
    if(level >= int(mipmaps.size())) level= int(mipmaps.size()) -1;
    
    assert(level < int(mipmaps.size()));
    const ImageData& image= mipmaps[level];
    return image.height;
}

int Texture::size( const Float lod ) const
{
    assert(mipmaps.size() > 0);
    
    int level= lod;
    if(level < base_level) level= base_level;
    if(level >= int(mipmaps.size())) level= int(mipmaps.size()) -1;
    
    assert(level < int(mipmaps.size()));
    const ImageData& image= mipmaps[level];
    return image.width * image.height;
}

Color Texture::sample( const Float u, const Float v, const Float lod ) const
{
    assert(mipmaps.size() > 0);
    
    int level= lod;
    if(level < base_level) level= base_level;
    if(level >= int(mipmaps.size())) level= int(mipmaps.size()) -1;
    
    assert(level >= 0);
    assert(level < int(mipmaps.size()));
    const ImageData& image= mipmaps[level];
    assert(image.width > 0);
    assert(image.height > 0);
    
#if 0    
    Float x= std::abs(std::fmod(std::floor(image.width * u), image.width));
    assert(x >= 0);
    assert(x < image.width);
    
    Float y= std::abs(std::fmod(std::floor(image.height * v), image.height));
    assert(y >= 0);
    assert(y < image.height);
#else
    
    int x= int(std::abs(std::floor(image.width * u))) % image.width;
    int y= int(std::abs(std::floor(image.height * v))) % image.height;
#endif

    size_t texel= image.offset(x, y);
    //~ assert(texel < image.data.size());
    
    Float alpha= 1;
    if(image.channels > 3) alpha= Float(image.data[texel +3]) / Float(255);
    
    return Color(Float(image.data[texel]) / Float(255), Float(image.data[texel +1]) / Float(255), Float(image.data[texel +2]) / Float(255), alpha);
}


TextureLib::TextureLib() : textures(), default_texture( read_image_data("data/debug2x2red.png") ), opaque(true) {}

int TextureLib::insert( const ImageData& data )
{
    assert(data.width > 0);
    assert(data.height > 0);
    textures.push_back(Texture(data));
    return (int) textures.size() -1;
}


TextureLib read_textures( std::vector<Material>& materials, const size_t max_size )
{
    TextureLib lib;
    
    // evalue la taille totale occuppee par toutes les images / textures
    size_t total_size= 0;
    for(int i= 0; i < int(materials.size()); i++)
    {
        Material& material= materials[i];
        
        if(!material.diffuse_filename.empty())
        {
            ImageData diffuse= read_image_data(material.diffuse_filename.c_str());
            if(diffuse.data.size())
            {
                total_size= total_size + miplevel_size(diffuse, 0);
                material.diffuse_texture= lib.insert(diffuse);
                
                // detecte les textures transparentes
                if(diffuse.channels == 4)
                {
                    const int w= diffuse.width;
                    const int h= diffuse.height;
                    const int stride= diffuse.channels * diffuse.size;
                    for(int y= 0; y < h; y++)
                    for(int x= 0; x < w; x++)
                        if(diffuse.data[offset(w, stride, x, y, 3)] < 255)
                        {
                            lib.opaque= false;
                            break;
                        }
                }
                
                // ne stocke plus les images apres avoir depasse la limite de taille
                if(total_size > max_size)
                    lib.textures[material.diffuse_texture].mipmaps[0].release();
            }
        }
        
        if(!material.specular_filename.empty())
        {
            ImageData specular= read_image_data(material.specular_filename.c_str());
            if(specular.data.size())
            {
                total_size= total_size + miplevel_size(specular, 0);
                material.specular_texture= lib.insert(specular);

                // ne stocke plus les images apres avoir depasse la limite de taille
                if(total_size > max_size)
                    lib.textures[material.specular_texture].mipmaps[0].release();
            }
        }
    }
    
    printf("%d textures.\n", lib.size());
    printf("  using %dMB / %dMB\n", int(total_size / 1024 / 1024), int(max_size / 1024 / 1024));
    if(lib.opaque)
        printf("  using opaque textures\n");
    else
        printf("  using alpha textures\n");
    
    // reduit les dimensions des images / textures jusqu'a respecter la limite de taille
    int lod= 0;
    for(; total_size > max_size ; lod++)
    {
        total_size= 0;
        for(int i= 0; i < lib.size(); i++)
            total_size= total_size + miplevel_size(lib.textures[i].mipmaps[0], lod);
    }
    printf("  lod %d, %dMB\n", lod, int(total_size / 1024 / 1024));
    
#if 1
    printf("building mipmaps...\n");
    // construit les textures a la bonne resolution
    for(int i= 0;  i < int(materials.size()); i++)
    {
        Material& material= materials[i];
        
        if(material.diffuse_texture != -1)
        {
            Texture& texture= lib.textures[material.diffuse_texture]; 
            if(texture.mipmaps[0].width > 0)
            {
                if(texture.mipmaps[0].data.empty())
                {
                    // recharge l'image, si necessaire
                    texture.mipmaps[0]= read_image_data(material.diffuse_filename.c_str());
                    assert(texture.mipmaps[0].data.size() > 0);
                }
                
                texture.build_mipmaps(lod);
                
                if(0)
                {
                    char tmp[1024];
                    for(int l= texture.base_level; l < int(texture.mipmaps.size()); l++)
                    {
                        sprintf(tmp, "%s%02d.png", material.diffuse_filename.c_str(), l);
                        write_image_data(texture.mipmaps[l], tmp);
                    }
                }
            }
        }
    
        if(materials[i].specular_texture != -1)
        {
            Texture& texture= lib.textures[material.specular_texture]; 
            if(texture.mipmaps[0].width > 0)
            {
                if(texture.mipmaps[0].data.empty())
                {
                    // recharge l'image, si necessaire
                    texture.mipmaps[0]= read_image_data(material.specular_filename.c_str());
                    assert(texture.mipmaps[0].data.size() > 0);
                }
                
                texture.build_mipmaps(lod);
            }
        }
    }

    if( 0 )
    {
        size_t mipmap_size= 0;
        for(int i= 0; i < int(materials.size()); i++)
        {
            if(materials[i].diffuse_texture != -1)
            {
                int id= materials[i].diffuse_texture;
                for(int k= 0; k < int(lib.textures[id].mipmaps.size()); k++)
                    mipmap_size+= lib.textures[id].mipmaps[k].data.size();
            }
            
            if(materials[i].specular_texture != -1)
            {
                int id= materials[i].specular_texture;
                for(int k= 0; k < int(lib.textures[id].mipmaps.size()); k++)
                    mipmap_size+= lib.textures[id].mipmaps[k].data.size();
            }
        }
        
        printf("  using %dMB / %dMB\n", int(mipmap_size / 1024 / 1024), int(max_size / 1024 / 1024));
    }
#endif

    return lib;
}


void release_textures( TextureLib& lib )
{
    return;
}

