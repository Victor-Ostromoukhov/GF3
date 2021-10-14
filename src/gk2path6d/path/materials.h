
#ifndef _MATERIALS_H
#define _MATERIALS_H

#include <vector>

#include "mesh.h"
#include "texturelib.h"


struct Materials
{
    TextureLib textures;
    std::vector<Material> materials;
    std::vector<unsigned int> triangles;
    
    Materials( ) {}
    
    const Material& material( const int triangle_id ) const
    {
        unsigned int id= triangles[triangle_id];
        return materials[id];
    }
    
    Color sample_diffuse( const int triangle_id, const Float u, const Float v, const Float lod ) const
    {
        unsigned int id= triangles[triangle_id];
        const Material& material= materials[id];
        if(material.diffuse_texture == -1)
            return White();
        
        assert(material.diffuse_texture < int(textures.size()));
        return textures.textures[material.diffuse_texture].sample(u, v, lod);
    }
    
    Color sample_specular( const int triangle_id, const Float u, const Float v, const Float lod ) const
    {
        unsigned int id= triangles[triangle_id];
        const Material& material= materials[id];
        if(material.specular_texture == -1)
            return White();
        
        assert(material.specular_texture < int(textures.size()));
        return textures.textures[material.specular_texture].sample(u, v, lod);
    }
};

int build_materials( Materials& materials, const Mesh& mesh );

#endif
