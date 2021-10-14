
#ifndef _STBVH_H
#define _STBVH_H

#include "bvh.h"


struct STBVHNode
{
    BVHBox bounds0;
    Float r0;
    Float r1;
    BVHBox bounds1;
    int left;
    int right;
    
    void make_node( const BVHBox& _b0, const Float _r0, const BVHBox& _b1, const Float _r1, const int _left, const int _right ) 
    { 
        bounds0= _b0; r0= _r0;
        bounds1= _b1; r1= _r1;
        left= _left;
        right= _right;
        assert(internal());
    }
    
    void make_leaf( const BVHBox& _b0, const Float _r0, const BVHBox& _b1, const Float _r1, const int _begin, const int _end ) 
    {
        bounds0= _b0; r0= _r0;
        bounds1= _b1; r1= _r1;
        left= _begin;
        right= - _end + _begin; // -n
        assert(leaf());
    }
    
    bool internal( ) const { return (right >= 0); }
    
    bool leaf( ) const { return (right < 0);}
    int leaf_begin( ) const { assert(leaf()); return left; }
    int leaf_end( ) const { assert(leaf()); return left - right; }
    int leaf_n( ) const { assert(leaf()); return - right; }
};



struct STBVHMesh 
{
    std::vector< const std::vector<vec3> * > vertices;
    std::vector<unsigned int> indices;
    int slots;
    int index;
    
    STBVHMesh( ) : vertices(), indices(), slots(0), index(-1) {}
    
    // time steps
    void buffers( const int n ) { slots= n; vertices.resize(n); }
    void shared_buffer( const int t, const std::vector<vec3>& data ) { assert(t < int(vertices.size())); vertices[t]= &data; }
    void index_buffer( const std::vector<unsigned int>& data ) { indices= data; }
    
    //
    int buffers( ) const { return slots; }
    const std::vector<vec3>& buffer( const int t ) const { assert(t < int(vertices.size())); return *vertices[t]; }
};


struct STBVHReference
{
    int object_id;
    int triangle_id;
    
    STBVHReference( const int _object, const int _triangle ) : object_id(_object), triangle_id(_triangle) {}
};

/* linear motion bvh,
    simplified stbvh 
    "STBVH: A Spatial-Temporal BVH for Efficient Multi-Segment Motion Blur", 
    woop hpg 2017
    http://www.sven-woop.de/papers/2017-HPG-msmblur.pdf

    slides http://www.highperformancegraphics.org/wp-content/uploads/2017/Papers-Session3/HPG2017_STBVH.pptx

    utilisation : creer des STBVHMesh puis les attacher aux STBVH, et construire le stbvh :
\code
    STBVH scene;
    
    // etape 1 : lire les keyframes
    std::vector<Mesh> frames;
    for(int i= 0; i < 50; i++)
        frames.push_back( read_mesh_cache( Format("bigguy/bigguy_%02d.obj", i)) );
    
    // etape 2 : ajouter les keyframes
    STBVHMesh bigguy;
    // nombre de keyframes
    bigguy.buffers(frames.size());
    
    // le vertex buffer de chaque key frame
    for(int i= 0; i < int(frames.size()); i++)
        // !! le vertex buffer doit rester accessible, le stbvh ne le copie pas, il le reference !!
        bigguy.shared_buffer(i, frames[i].positions());
    
    // l'index buffer (commun a toutes les keyframes)
    bigguy.index_buffer(frames[0].indices());
    // inserer l'objet dans le stbvh
    scene.attach(bigguy);
    
    // etape 3 : ajouter un objet statique
    Mesh ground_mesh= read_mesh_cache( "bigguy/ground.obj" );

    STBVHMesh ground;
    ground.buffers(1);
    ground.shared_buffer(0, ground_mesh.positions());
    ground.index_buffer(ground_mesh.indices());
    scene.attach(ground);
    
    // etape 4 : construire le stbvh
    STBVHBuilder tmp;
    tmp.build(scene);
\endcode
 */

struct IntersectionFilter
{
    IntersectionFilter( ) {}
    virtual ~IntersectionFilter( ) {}
    
    virtual bool operator() ( const Ray& ray, const Hit& hit ) const {  return true; }
};


struct STBVH
{
    std::vector<STBVHMesh> objects;
    std::vector<STBVHReference> references;
    std::vector<STBVHNode> nodes;
    BVHBox bounds;
    int root;
    bool notime;
    
    STBVH( ) : objects(), references(), nodes(), bounds(), root(-1), notime(false) {}
    
    int attach( const STBVHMesh& object ) 
    { 
        assert(object.index == -1); 
        assert(object.slots > 0); 
        assert(object.indices.size()); 
        
        int index= int(objects.size()); 
        objects.push_back(object); 
        objects.back().index= index; 
        return index;
    }
    
    //
    bool intersect( const Ray& ray, Hit& hit, const IntersectionFilter& filter ) const
    {
        Vector invd(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);

        hit.t= ray.tmax;
        hit.triangle_id= -1;
        hit.mesh_id= -1;
        
        if(notime) intersect_notime(root, ray, invd, hit, filter);
        else intersect(root, ray, invd, hit);
        return (hit.mesh_id != -1);
    }

    bool visible( const Float time, const Point& p, const Point& q, const IntersectionFilter& filter ) const
    {
        Ray ray(time, p, q);
        Vector invd(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        
        if(notime) return visible_notime(root, ray, invd, filter);
        else return visible(root, ray, invd);
    }
    
    bool visible( const Float time, const Point& p, const Vector& d, const IntersectionFilter& filter ) const
    {
        Ray ray(time, p, d);
        Vector invd(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        
        if(notime) return visible_notime(root, ray, invd, filter);
        else return visible(root, ray, invd);
    }
    
protected:
    void intersect( const int index, const Ray& ray, const Vector& invd, Hit& hit ) const;
    bool visible( const int index, const Ray& ray, const Vector& invd ) const;
    
    void intersect_notime( const int index, const Ray& ray, const Vector& invd, Hit& hit, const IntersectionFilter& filter ) const;
    bool visible_notime( const int index, const Ray& ray, const Vector& invd, const IntersectionFilter& filter ) const;
};


#endif
