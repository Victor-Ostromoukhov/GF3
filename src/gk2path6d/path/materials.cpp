
#include "materials.h"


int build_materials( Materials& materials, const Mesh& mesh )
{
    materials.materials= mesh.mesh_materials();
    materials.triangles= mesh.materials();
    printf("%d materials.\n", (int) materials.materials.size());
    
    materials.textures= read_textures(materials.materials);
    return (int) materials.materials.size();
}
