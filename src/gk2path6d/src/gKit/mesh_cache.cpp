// stat
#ifdef _MSC_VER
    #include <sys/types.h>
#endif
#include <sys/stat.h>

#include <cassert>

#include "mesh.h"
#include "wavefront.h"
#include "mesh_buffer.h"
#include "mesh_cache.h"

bool infos( const char *filename, size_t& time )
{
    time= 0;
    
#ifndef _MSC_VER
    struct stat info;
    if(stat(filename, &info) < 0)
        return false;

    //~ size= info.st_size;
    time= info.st_mtime;
    return true;

#else
    struct _stat64 info;
    if(_stat64(filename, &info) < 0)
        return false;

    //~ size= (size_t) info.st_size;
    time = (size_t) info.st_mtime;
    return true;
#endif
}

template< typename T>
static int read( FILE *in, std::vector<T>& data )
{
    data.clear();
    
    int n= 0;
    size_t code= fread(&n, sizeof(int), 1, in);
    assert(code == 1);
    if(code == 0 || n == 0)
        return 0;
    
    data.resize(n);
    code= fread(data.data(), sizeof(T), data.size(), in);
    assert(code == n);
    return n;
}

Mesh read_mesh_cache( const char *filename )
{
    Mesh mesh;
    size_t obj_time;
    if(!infos(filename, obj_time))
    {
        printf("[error] loading mesh '%s'...\n", filename);
        return mesh;
    }
    
    size_t cache_time;
    std::string cache_filename= std::string(filename) + ".mesh";
    if(infos(cache_filename.c_str(), cache_time) && obj_time < cache_time)
    {
        FILE *in= fopen(cache_filename.c_str(), "rb");
        assert(in);

        printf("loading mesh cache '%s'...\n", cache_filename.c_str());
        
        int primitives= mesh.primitives();
        int code = fread(&primitives, sizeof(int), 1, in);
        assert(code == 1);
        
        mesh= Mesh(primitives);
        
        std::vector<vec3> data3;
        if(read(in, data3)) 
            mesh.positions(data3);
        
        std::vector<vec2> data2;
        if(read(in, data2))
            mesh.texcoords(data2);
        if(read(in, data3))
            mesh.normals(data3);
        
        std::vector<vec4> data4;
        if(read(in, data4))
            mesh.colors(data4);
        
        std::vector<unsigned int> data;
        if(read(in, data))
            mesh.indices(data);
        if(read(in, data))
            mesh.materials(data);
        fclose(in);
        
        cache_filename= std::string(filename) + ".materials";
        MaterialLib materials= read_materials(cache_filename.c_str());
        mesh.mesh_materials(materials.data);
        
        printf("  mesh: %d positions, %d texcoords, %d normals, %d indices, %d materials\n", 
            (int) mesh.positions().size(), (int) mesh.texcoords().size(), (int) mesh.normals().size(), (int) mesh.indices().size(), (int) mesh.mesh_materials().size());
        
        //~ printf("    indices %dMB\n", int(mesh.indices().size() * sizeof(int) / 1024 / 1024));
        //~ printf("    positions %dMB\n", int(mesh.positions().size() * sizeof(vec3) / 1024 / 1024));
        //~ printf("          %dMB\n", int((mesh.indices().size() * sizeof(int) + mesh.positions().size() * sizeof(vec3)) / 1024 / 1024));
        //~ printf("    array %dMB\n", 3 * int(mesh.indices().size() * sizeof(vec3) / 1024 / 1024));
    }
    else
    {
        MeshData data= read_mesh_data(filename);
        mesh= build_linear_mesh(data);
        
        if(!mesh.triangle_count())
            return mesh;
        
        write_mesh_cache(mesh, filename);
    }
    
    return mesh;
}


template< typename T >
static size_t write( const std::vector<T>& buffer, FILE *out )
{
    int n= int(buffer.size());
    size_t code= fwrite(&n, sizeof(int), 1, out);
    if(n == 0)
        return 1;
    assert(code == 1);

    code= fwrite(buffer.data(), sizeof(T), buffer.size(), out);
    assert(code == buffer.size());
    return code;
}

int write_mesh_cache( const Mesh& mesh, const char *filename )
{
    if(!mesh.triangle_count())
        return -1;
    
    std::string cache_filename= std::string(filename) + ".mesh";
    FILE *out= fopen(cache_filename.c_str(), "wb");
    if(!out)
        return -1;
    
    // mesh data
    printf("writing mesh cache '%s'...\n", cache_filename.c_str());
    
    int primitives= mesh.primitives();
    int code= fwrite(&primitives, sizeof(int), 1, out);
    assert(code == 1);
    
    write(mesh.positions(), out);
    write(mesh.texcoords(), out);
    write(mesh.normals(), out);
    write(mesh.colors(), out);
    write(mesh.indices(), out);
    
    write(mesh.materials(), out);
    code= fclose(out);
    assert(code == 0);
    
    // materials
    cache_filename= std::string(filename) + ".materials";
    out= fopen(cache_filename.c_str(), "wt");
    if(!out)
        return -1;
    
    printf("writing materials '%s'...\n", cache_filename.c_str());
    
    const std::vector<Material>& materials= mesh.mesh_materials();
    for(int i= 0; i < int(materials.size()); i++)
    {
        const Material& m= materials[i];
        
        fprintf(out, "newmtl %d\n", i);
        if(m.diffuse.power())
            fprintf(out, "Kd %lf %lf %lf\n", double(m.diffuse.r), double(m.diffuse.g), double(m.diffuse.b));
        if(m.specular.power())
            fprintf(out, "Ks %lf %lf %lf\n", double(m.specular.r), double(m.specular.g), double(m.specular.b));
        if(m.emission.power())
            fprintf(out, "Ke %lf %lf %lf\n", double(m.emission.r), double(m.emission.g), double(m.emission.b));
        if(m.ns)
            fprintf(out, "Ns %lf\n", double(m.ns));
        if(!m.diffuse_filename.empty())
            fprintf(out, "map_Kd %s\n", m.diffuse_filename.c_str());
        if(!m.specular_filename.empty())
            fprintf(out, "map_Ks %s\n", m.specular_filename.c_str());
        fprintf(out, "\n");
    }
    fclose(out);
    
    return 0;
}

