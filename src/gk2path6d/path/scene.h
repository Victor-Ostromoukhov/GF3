
#ifndef _SCENE_H
#define _SCENE_H

#include <vector>
#include <cassert>
#include <cmath>

#include "vec.h"
#include "color.h"
#include "ray.h"
#include "triangle.h"

#include "mesh.h"
#include "mesh_cache.h"

#include "stbvh.h"
#include "stbvh_builder.h"

#include "materials.h"
#include "points.h"


// transforme les positions / normales des triangles du mesh.
Mesh transform( const Mesh& mesh, const Transform& m );


struct SceneObject
{
    Materials materials;
    std::vector<Mesh> frames;
    bool noalpha;
    
    SceneObject( ) : materials(), frames(), noalpha(true) {}
    SceneObject( const char *prefix, const int start, const int stop, const Transform& m= Identity() ) : materials(), frames(), noalpha(true) { read_frames(prefix, start, stop, m); }
    SceneObject( const char *filename, const Transform& m= Identity() ) : materials(), frames(), noalpha(true) { read_mesh(filename, m); }
    
    void read_frames( const char *prefix, const int start, const int stop, const Transform& m= Identity() )
    {
        assert(frames.size() == 0);
        
        char tmp[4096];
        for(int i= start; i <= stop; i++)
        {
            sprintf(tmp, prefix, i);
            frames.push_back( transform(read_mesh_cache(tmp), m) );
        }
    }
    
    void read_mesh( const char *filename, const Transform& m= Identity() )
    {
        frames.push_back( transform(read_mesh_cache(filename), m) );
    }

    // recuperer les sources de lumiere 
    // todo sources animees
    void build_sources( std::vector<Source>& sources )
    {
        for(int i= 0; i < frames[0].triangle_count(); i++)
        {
            // recupere la matiere associee a chaque triangle de l'objet
            Color emission= frames[0].triangle_material(i).emission;
            Float area= Triangle(frames[0].triangle(i)).area();
            
            if(emission.power() * area > 0)
                sources.push_back( Source(frames[0].triangle(i), emission) );
        }
    }
    
    int frame0( const Float time ) const
    { 
        int f0= std::floor(time * Float(frames.size() - 1)); 
        assert(f0 >= 0); 
        assert(f0 < int(frames.size())); 
        return f0;
    }
    
    int frame1( const Float time ) const
    { 
        int f1= std::ceil(time * Float(frames.size() - 1)); 
        assert(f1 >= 0); 
        assert(f1 < int(frames.size())); 
        return f1;
    }
    
    vec2 texcoords( const Raydata& ray ) const
    {
        int t0= frame0(ray.time());
        int t1= frame1(ray.time());
        Float t= ray.time() * Float(frames.size() -1) - t0;
        if(t < 0) t= 0;
        if(t > 1) t= 1;
        
        Triangle triangle0= Triangle(frames[t0].triangle(ray.hitid()));
        Triangle triangle1= Triangle(frames[t1].triangle(ray.hitid()));
        
        vec2 texcoord0= triangle0.texcoord(ray.hitu(), ray.hitv());
        vec2 texcoord1= triangle1.texcoord(ray.hitu(), ray.hitv());
        return vec2(texcoord0.x * (1 - t) + texcoord1.x * t, texcoord0.y * (1 - t) + texcoord1.y * t);
    }
    
    Vector normal( const Raydata& ray ) const
    {
        int t0= frame0(ray.time());
        int t1= frame1(ray.time());
        Float t= ray.time() * Float(frames.size() -1) - t0;
        if(t < 0) t= 0;
        if(t > 1) t= 1;
        
        Triangle triangle0= Triangle(frames[t0].triangle(ray.hitid()));
        Triangle triangle1= Triangle(frames[t1].triangle(ray.hitid()));
        
        Vector normal0= triangle0.normal(ray.hitu(), ray.hitv());
        Vector normal1= triangle1.normal(ray.hitu(), ray.hitv());
        return normalize(normal0 * (1 - t) + normal1 * t);
    }
    
    Color diffuse_texture( const Raydata& ray, const Float lod= 0 ) const
    { 
        vec2 uv= texcoords(ray);
        return materials.sample_diffuse(ray.hitid(), uv.x, uv.y, lod);
    }
    
    Color specular_texture( const Raydata& ray, const Float lod= 0 ) const
    {
        vec2 uv= texcoords(ray);
        return materials.sample_specular(ray.hitid(), uv.x, uv.y, lod);
    }
    
    const Material& material( const Raydata& ray ) const
    {
        return materials.material(ray.hitid());
    }
};


struct Scene
{
    STBVH bvh;
    
    std::vector<SceneObject> objects;
    int max_frames;
    
    std::string camera_filename;
    std::string lens_filename;
    
    Scene( ) : bvh(), objects(), max_frames(0), camera_filename(), lens_filename() {}
    
    int attach( const SceneObject& object ) 
    {
        assert(object.frames.size());
        max_frames= std::max(max_frames, int(object.frames.size()));
        
        objects.push_back(object);
        return int(objects.size()) -1;
    }
    
    void build( )
    {
        assert(objects.size());
        
        for(int i= 0; i < int(objects.size()); i++)
        {
            build_materials(objects[i].materials, objects[i].frames[0]);
            
            STBVHMesh mesh;
            const SceneObject& object= objects[i];
            
            mesh.buffers(object.frames.size());
            for(int i= 0; i < int(object.frames.size()); i++)
                mesh.shared_buffer(i, object.frames[i].positions());
            
            mesh.index_buffer(object.frames[0].indices());
            bvh.attach(mesh);
            // plus besoin de mesh
        }
        
        printf("building STBVH...\n");
        STBVHBuilder tmp;
        tmp.build(bvh);
    }
    
    void build_sources( std::vector<Source>& sources )
    {
        for(int i= 0; i < int(objects.size()); i++)
            objects[i].build_sources(sources);
    }
    
    bool intersect( Raydata& data ) const;
    //~ {
        //~ IntersectionAlpha filter(*this, 0.5);
        //~ bvh.intersect(data.ray, data.hit, filter);
        //~ return data.intersect();
    //~ }
    
    bool visible( const Float t, const Point& o, const Point& e ) const;
    //~ {
        //~ return bvh.visible(t, o, e);
    //~ }
    
    bool visible( const Float t, const Point& o, const Vector& d ) const;
    //~ {
        //~ return bvh.visible(t, o, d);
    //~ }    
    
    int frames( ) const { return max_frames; }
    int frame0( const Float time ) const { return std::floor(time * Float(max_frames - 1));  }
    int frame1( const Float time ) const { return std::ceil(time * Float(max_frames - 1));  }
    
    vec2 texcoords( const Raydata& ray ) const
    {
        assert(ray.hit.mesh_id != -1);
        assert(ray.hit.mesh_id < int(objects.size()));
        return objects[ray.hit.mesh_id].texcoords(ray);
    }
    
    Vector normal( const Raydata& ray ) const
    {
        assert(ray.hit.mesh_id != -1);
        assert(ray.hit.mesh_id < int(objects.size()));
        return objects[ray.hit.mesh_id].normal(ray);
    }
    
    const Material& material( const Raydata& ray ) const
    {
        assert(ray.hit.mesh_id != -1);
        assert(ray.hit.mesh_id < int(objects.size()));
        return objects[ray.hit.mesh_id].material(ray);
    }
    
    Color diffuse_texture( const Raydata& ray, const Float lod ) const
    { 
        assert(ray.hit.mesh_id != -1);
        assert(ray.hit.mesh_id < int(objects.size()));
        return objects[ray.hit.mesh_id].diffuse_texture(ray, lod);
    }
    
    Color specular_texture( const Raydata& ray, const Float lod= 0 ) const    
    {
        assert(ray.hit.mesh_id != -1);
        assert(ray.hit.mesh_id < int(objects.size()));
        return objects[ray.hit.mesh_id].specular_texture(ray, lod);
    }
};


struct IntersectionAll : public IntersectionFilter
{
    IntersectionAll( ) : IntersectionFilter() {}
    
    bool operator() ( const Ray& ray, const Hit& hit ) const
    {
        return true;
    }
};

struct IntersectionAlpha : public IntersectionFilter
{
    const Scene& scene;
    Float alpha;
    
    IntersectionAlpha( const Scene& s, const Float a= Float(0.3) ) : IntersectionFilter(), scene(s), alpha(a) {}
    
    bool operator() ( const Ray& ray, const Hit& hit ) const
    {
        Raydata data(ray, hit);
        Color color= scene.diffuse_texture(data, 0);
        return (color.a > alpha);       // surface opaque, accepter l'intersection
    }
};


Scene read_scene( const char *filename );

#endif
