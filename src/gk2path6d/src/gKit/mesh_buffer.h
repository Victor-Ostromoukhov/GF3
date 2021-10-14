
#ifndef _MESH_BUFFER_H
#define _MESH_BUFFER_H

#include "mesh.h"
#include "wavefront.h"


//! representation d'une sequence de triangles associes a la meme matiere
struct MeshGroup
{
    int material;       //!< indice de la matiere
    int first;          //!< indice des premiers sommets
    int count;          //!< nombre d'indices
    
    MeshGroup( const int _id= -1, const int _first= 0 ) : material(_id), first(_first), count(0) {}
};


//! representation d'un objet.
struct MeshBuffer
{
    std::vector<vec3> positions;                //!< attribut position
    std::vector<vec2> texcoords;                //!< attribut coordonnees de texture
    std::vector<vec3> normals;                  //!< attribut normale
    
    std::vector<int> material_indices;          //!< indice de la matiere des triangles
    std::vector<int> indices;                   //!< indices des sommets des triangles
    
    std::vector<Material> materials;        //!< ensemble de matieres
    std::vector<MeshGroup> material_groups;     //!< sequence de triangles groupes par matiere
};


//! construction a partir des donnees d'un fichier .obj.
MeshBuffer build_groups( const MeshData& data );
Mesh build_linear_mesh( const MeshData& data );

void build_unique( MeshData& data );

//! (re-) calcule les normales des sommets. utiliser avant les reindexations, cf indices() et vertices().
void build_normals( MeshData& data );

//! construit les triangles. prepare l'affichage openGL, avec glDrawArrays().
Mesh build_triangle_mesh( const MeshData& data );

#endif
