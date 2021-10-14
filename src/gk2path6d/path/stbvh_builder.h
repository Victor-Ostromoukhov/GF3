
#ifndef _STBVH_BUILDER_H
#define _STBVH_BUILDER_H

#include "stbvh.h"


struct STBVHPrimitive
{
    BVHBox bounds0;     // linear bounds
    BVHBox bounds1;
    
    BVHBox bounds05;
    Point center05;
    
    int object;         // bvh.objects[object]
    int index;          // triangle(indices[index], indices[index +1], indices[index +2])
    
    int n;              // motion segments
    int t0;             // first time slot : object.vertices[t0]
    int t1;             // last time slot : object.vertices[t1]
    
    STBVHPrimitive( ) : object(-1), index(-1), n(0), t0(0), t1(0) {}
    STBVHPrimitive( const BVHBox& b0, const BVHBox& b1, const int _object, const int _id, const int _n, const int _t0, const int _t1 )
        : 
        bounds0(b0), bounds1(b1), bounds05(bounds(.5f)), center05(bounds05.center()), 
        object(_object), index(_id), n(_n), t0(_t0), t1(_t1)
    {
        assert(n > 1);
        assert(t0 <= t1);
    }
    
    // interpole la bbox pour l'instant t
    BVHBox bounds( const Float t ) const { return BVHBox( bounds0.pmin * (1 - t) + bounds1.pmin * t, bounds0.pmax * (1 - t) + bounds1.pmax * t ); }
    
    // renvoie le nombre de segments inclus dans l'intervalle r0 r1
    int motion( const Float r0, const Float r1 ) { return slot1(r1) - slot0(r0); }
    
    // renvoie l'indice d'un segment pour l'instant r0 < r
    int slot0( const Float r ) { return std::floor(r * Float(n -1)); }
    
    // renvoie l'indice d'un segment pour l'instant r1 >= r
    int slot1( const Float r ) { return std::ceil(r * Float(n -1)); }
    
    Float round0( const Float r ) { return slot0(r) / Float(n -1); }
    Float round1( const Float r ) { return slot1(r) / Float(n -1); }
};


/* linear motion bvh,
    simplified stbvh 
    "STBVH: A Spatial-Temporal BVH for Efficient Multi-Segment Motion Blur", 
    woop hpg 2017
    http://www.sven-woop.de/papers/2017-HPG-msmblur.pdf

    slides http://www.highperformancegraphics.org/wp-content/uploads/2017/Papers-Session3/HPG2017_STBVH.pptx
 */

struct STBVHBuilder
{
    void build( STBVH& bvh );

protected:
    std::vector<STBVHNode> nodes;
    std::vector<STBVHReference> references;
    
    //! internal build procedure.
    int build_node( const STBVH& bvh, std::vector<STBVHPrimitive>& primitives, const Float r0, const Float r1 );

    STBVHPrimitive build_primitive( const STBVHMesh& object, const int index, const int t0, const int t1 );
    void build_bounds( const STBVH& bvh, std::vector<STBVHPrimitive>& primitives, const Float r0, const Float r1, BVHBox& b0, BVHBox& b1 );
};

#endif
