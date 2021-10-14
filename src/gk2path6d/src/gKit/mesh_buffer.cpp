
#include <cassert>
#include <map>
#include <algorithm>

#include "mesh_buffer.h"


// compare la matiere de 2 triangles
struct compareMaterial
{
    const std::vector<unsigned int>& material_buffer;
    
    compareMaterial( const std::vector<unsigned int>& _buffer ) : material_buffer(_buffer) {}
    
    bool operator() ( const int& a, const int& b ) const
    {
        return material_buffer[a] < material_buffer[b];
    }
};


//! representation de l'indexation complete d'un sommet
struct MeshVertex
{
    int material;
    int position;
    int texcoord;
    int normal;
    
    MeshVertex( ) : material(-1), position(-1), texcoord(-1), normal(-1) {}
    MeshVertex( const int m, const int p, const int t, const int n ) : material(m), position(p), texcoord(t), normal(n) {}
    
    // comparaison lexicographique de 2 indices
    bool operator< ( const MeshVertex& b ) const
    {
        if(material != b.material) return material < b.material;
        if(position != b.position) return position < b.position;
        if(texcoord != b.texcoord) return texcoord < b.texcoord;
        if(normal != b.normal) return normal < b.normal;
        return false;
    }
};


MeshBuffer build_groups( const MeshData& data )
{
    MeshBuffer mesh;

    mesh.materials= data.materials;
    mesh.material_indices.reserve(data.material_indices.size());
    
    mesh.positions.reserve(data.positions.size());
    mesh.texcoords.reserve(data.positions.size());
    mesh.normals.reserve(data.positions.size());
    
    // tri les triangles par matiere
    std::vector<int> triangles;
    triangles.reserve(data.material_indices.size());
    for(int i= 0; i < (int) data.material_indices.size(); i++)
        triangles.push_back(i);
    
    std::stable_sort(triangles.begin(), triangles.end(), compareMaterial(data.material_indices));
    
    // groupes de triangles 
    unsigned int material_id= data.material_indices[triangles[0]];
    mesh.material_groups.push_back( MeshGroup(material_id, 0) );
    
    bool has_texcoords= !data.texcoords.empty();
    bool has_normals= !data.normals.empty();
    
    // re-ordonne les triangles et les attributs
    std::map<MeshVertex, int> remap;
    for(int i= 0; i < (int) triangles.size(); i++)
    {
        // matiere du triangle
        mesh.material_indices.push_back(data.material_indices[triangles[i]]);
        
        // associe la matiere de la face a ses sommets
        if(material_id != data.material_indices[triangles[i]])
        {
            // construit les groupes de triangles associes a la meme matiere
            material_id= data.material_indices[triangles[i]];
            
            mesh.material_groups.back().count= 3*i - mesh.material_groups.back().first;
            mesh.material_groups.push_back( MeshGroup(material_id, 3*i) );
        }
        
        for(int k= 0; k < 3; k++)
        {
            // indice du kieme sommet du ieme triangle re-ordonne
            int index= 3*triangles[i] + k;
            MeshVertex vertex= MeshVertex(material_id, data.position_indices[index], has_texcoords ? data.texcoord_indices[index] : -1, has_normals ? data.normal_indices[index] : -1);
            
            auto found= remap.insert( std::make_pair(vertex, int(remap.size())) );
            if(found.second)
            {
                // copie les attributs
                assert(data.position_indices[index] != -1);
                mesh.positions.push_back( data.positions[data.position_indices[index]] );
                
                if(data.texcoord_indices[index] != -1)
                    // copie les texcoord du sommet, si elles sont definies
                    mesh.texcoords.push_back( data.texcoords[data.texcoord_indices[index]] );
                else if(has_texcoords)
                    // copie une valeur par defaut, tous les sommets n'ont pas de texcoord
                    mesh.texcoords.push_back( vec2() );
                
                if(data.normal_indices[index] != -1) 
                    // copie la normale du sommet, si ell est definie
                    mesh.normals.push_back( data.normals[data.normal_indices[index]] );
                else if(has_normals)
                    // copie une valeur par defaut, tous les sommets n'ont pas de normale
                    mesh.normals.push_back( vec3() );
            }
            
            // construit l'index buffer
            mesh.indices.push_back(found.first->second);
        }
    }

    // termine la description du dernier groupe de triangles
    mesh.material_groups.back().count= 3*triangles.size() - mesh.material_groups.back().first;
    
    printf("buffers: %d positions, %d texcoords, %d normals, %d indices, %d groups\n", 
        (int) mesh.positions.size(), (int) mesh.texcoords.size(), (int) mesh.normals.size(), (int) mesh.indices.size(), (int) mesh.material_groups.size());
    
    return mesh;
}


Mesh build_linear_mesh( const MeshData& data )
{
    // attributs
    std::vector<vec3> positions;
    positions.reserve(data.positions.size());
    std::vector<vec2> texcoords;
    texcoords.reserve(data.positions.size());
    std::vector<vec3> normals;
    normals.reserve(data.positions.size());
    
    // matieres
    std::vector<unsigned int> indices;
    indices.reserve(data.material_indices.size());
    std::vector<unsigned int> material_indices;
    material_indices.reserve(data.material_indices.size());
    
    // tri les triangles par matiere
    std::vector<int> triangles;
    triangles.reserve(data.material_indices.size());
    for(int i= 0; i < (int) data.material_indices.size(); i++)
        triangles.push_back(i);
    
    std::stable_sort(triangles.begin(), triangles.end(), compareMaterial(data.material_indices));
    
    bool has_texcoords= !data.texcoords.empty();
    bool has_normals= !data.normals.empty();
    
    // re-ordonne les triangles et les attributs
    std::map<MeshVertex, int> remap;
    for(int i= 0; i < int(triangles.size()); i++)
    {
        // matiere du triangle
        int material_id=  data.material_indices[triangles[i]];
        material_indices.push_back(material_id);
        
        for(int k= 0; k < 3; k++)
        {
            // indice du kieme sommet du ieme triangle re-ordonne
            int index= 3*triangles[i] + k;
            //~ MeshVertex vertex= MeshVertex(material_id, data.position_indices[index], data.texcoord_indices[index], data.normal_indices[index]);
            MeshVertex vertex= MeshVertex(-1, data.position_indices[index], has_texcoords ? data.texcoord_indices[index] : -1, has_normals ? data.normal_indices[index] : -1);
            
            auto found= remap.insert( std::make_pair(vertex, int(remap.size())) );
            if(found.second)
            {
                // copie les attributs
                assert(data.position_indices[index] != -1);
                positions.push_back( data.positions[data.position_indices[index]] );
                
                if(has_texcoords && data.texcoord_indices[index] != -1)
                    // copie les texcoord du sommet, si elles sont definies
                    texcoords.push_back( data.texcoords[data.texcoord_indices[index]] );
                else if(has_texcoords)
                    // copie une valeur par defaut, tous les sommets n'ont pas de texcoord
                    texcoords.push_back( vec2() );
                
                if(has_normals && data.normal_indices[index] != -1) 
                    // copie la normale du sommet, si ell est definie
                    normals.push_back( data.normals[data.normal_indices[index]] );
                else if(has_normals)
                    // copie une valeur par defaut, tous les sommets n'ont pas de normale
                    normals.push_back( vec3() );
            }
            
            // construit l'index buffer
            indices.push_back(found.first->second);
        }
    }

    Mesh mesh(GL_TRIANGLES);    // \todo vire la dependance openGL
    mesh.mesh_materials(data.materials);
    mesh.materials(material_indices);
    
    mesh.positions(positions);
    mesh.texcoords(texcoords);
    mesh.normals(normals);
    mesh.indices(indices);
    
    printf("  linear mesh: %d positions, %d texcoords, %d normals, %d indices, %d materials\n", 
        (int) mesh.positions().size(), (int) mesh.texcoords().size(), (int) mesh.normals().size(), (int) mesh.indices().size(), (int) mesh.mesh_materials().size());
    
    return mesh;
}


void bounds( const MeshData& data, Point& pmin, Point& pmax )
{
    if(data.positions.size() < 1)
        return;
    
    pmin= Point(data.positions[0]);
    pmax= pmin;

    for(int i= 1; i < int(data.positions.size()); i++)
    {
        vec3 p= data.positions[i];
        pmin= Point( std::min(pmin.x, p.x), std::min(pmin.y, p.y), std::min(pmin.z, p.z) );
        pmax= Point( std::max(pmax.x, p.x), std::max(pmax.y, p.y), std::max(pmax.z, p.z) );
    }
}


struct vec3_less
{
    bool operator( ) ( const vec3& a, const vec3& b ) const
    {
        if(a.x != b.x) return a.x < b.x;
        if(a.y != b.y) return a.y < b.y;
        if(a.z != b.z) return a.z < b.z;
        
        return false;
    }
};

void build_unique( MeshData& data )
{
    std::vector<vec3> positions;
    std::map<vec3, int, vec3_less> remap;
    for(int i= 0; i < int(data.positions.size()); i++)
    {
        auto found= remap.insert( std::make_pair(data.positions[i], int(positions.size())) );
        if(found.second)
            positions.push_back(data.positions[i]);
    }
    printf("  positions %d, unique %d\n", int(data.positions.size()), int(remap.size()));
    
    for(int i= 0; i < int(data.position_indices.size()); i++)
    {
        vec3 p= data.positions[data.position_indices[i]];
        data.position_indices[i]= remap[p];
    }
    
    std::swap(data.positions, positions);
}


void build_normals( MeshData& data )
{
    // une normale par position
    std::vector<Vector> normals(data.positions.size(), Vector());
    
    for(int i= 0; i + 2 < int(data.position_indices.size()); i+= 3)
    {
        // positions des sommets du triangle
        int v0= data.position_indices[i];
        int v1= data.position_indices[i +1];
        int v2= data.position_indices[i +2];
        
        Point a= Point(data.positions[v0]);
        Point b= Point(data.positions[v1]);
        Point c= Point(data.positions[v2]);
        
        // normale geometrique
        Vector n= normalize(cross(normalize(b - a), normalize(c - a)));
        
        // somme la normale sur les sommets du triangle
        normals[v0]= normals[v0] + n;
        normals[v1]= normals[v1] + n;
        normals[v2]= normals[v2] + n;
    }
    
    // copie / conversion
    data.normals.clear();
    data.normals.reserve(normals.size());
    for(int i= 0; i < int(normals.size()); i++)
        data.normals.push_back( vec3(normalize(normals[i])) );
    
    // re-indexe les sommets
    for(int i= 0; i < int(data.normal_indices.size()); i++)
        data.normal_indices[i]= data.position_indices[i];
}


Mesh build_triangle_mesh( const MeshData& data )
{
    std::vector<vec3> positions;
    positions.reserve(data.positions.size());
    std::vector<vec2> texcoords;
    texcoords.reserve(data.texcoords.size());
    std::vector<vec3> normals;
    normals.reserve(data.normals.size());
    
    bool has_normals= (data.position_indices.size() == data.normal_indices.size());
    bool has_texcoords= (data.position_indices.size() == data.texcoord_indices.size());
    for(int i= 0; i < int(data.position_indices.size()); i++)
    {
        positions.push_back( data.positions[data.position_indices[i]] );
        
        if(has_texcoords && data.texcoord_indices[i] != -1)
            texcoords.push_back( data.texcoords[data.texcoord_indices[i]] );
        else if(has_texcoords)
            texcoords.push_back( vec2() );
        
        if(has_normals && data.normal_indices[i] != -1)
            normals.push_back( data.normals[data.normal_indices[i]] );
        else if(has_normals)
            normals.push_back( vec3() );
    }
    
    Mesh mesh(GL_TRIANGLES);    // \todo vire la dependance openGL
    mesh.mesh_materials(data.materials);
    mesh.materials(data.material_indices);
    mesh.positions(positions);
    mesh.texcoords(texcoords);
    mesh.normals(normals);
    
    return mesh;
}
