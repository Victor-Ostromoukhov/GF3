
#ifndef _OBJ_H
#define _OBJ_H

#include "mesh.h"


//! \addtogroup objet3D
///@{

//! \file 
//! charge un fichier wavefront .obj et construit un mesh.

//! charge un fichier wavefront .obj et renvoie un mesh compose de triangles non indexes. utiliser glDrawArrays pour l'afficher. a detruire avec Mesh::release( ).
Mesh read_mesh( const char *filename );

//! enregistre un mesh dans un fichier .obj.
int write_mesh( const Mesh& mesh, const char *filename );

//! ensemble de matieres.
struct MaterialLib
{
    std::vector<std::string> names;
    std::vector<Material> data;
};

//! charge une description de matieres, utilise par read_mesh.
MaterialLib read_materials( const char *filename );


//! stockage temporaire des donnees d'un ficher .obj.
struct MeshData
{
    std::vector<vec3> positions;
    std::vector<vec2> texcoords;
    std::vector<vec3> normals;
    
    std::vector<int> position_indices;
    std::vector<int> texcoord_indices;
    std::vector<int> normal_indices;
    
    std::vector<Material> materials;
    std::vector<unsigned int> material_indices;
};


//! charge un fichier wavefront .obj et renvoie les donnees.
MeshData read_mesh_data( const char *filename );

int write_mesh_data( const MeshData& data, const char *filename );

//! utilitaires manipulation de noms de fichiers
std::string pathname( const std::string& filename );
std::string normalize_path( const std::string& filename );
bool relative_path( const std::string& filename );

///@}
#endif
