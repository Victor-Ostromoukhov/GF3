
#ifndef _MESH_CACHE_H
#define _MESH_CACHE_H

#include "mesh.h"

bool infos( const char *filename, size_t& time );

Mesh read_mesh_cache( const char *filename );
int write_mesh_cache( const Mesh& mesh, const char *filename );

#endif
