
#ifndef _BVH_BUILDER_H
#define _BVH_BUILDER_H

#include "bvh.h"


//! functor interface: evaluates the cost a of cut.
struct BVHCost
{
    Float node_cost;
    Float triangle_cost;
    
    BVHCost( const Float _node= 1, const Float _triangle= 1 ) :  node_cost(_node), triangle_cost(_triangle) {}
    virtual ~BVHCost( ) {}
    
    //! cost of subdividing a node.
    virtual Float operator() ( const BVHBox& bbox_left, const int triangle_left, const BVHBox& bbox_right, const int triangle_right ) const = 0;
    
    //! leaf cost / cost of not subdividing a node.
    virtual Float operator() ( const int triangles ) const = 0;
};


//! functor: evaluates the SAH cost.
struct SAHCost : public BVHCost
{
    SAHCost( const Float _node= 1, const Float _triangle= 1 ) : BVHCost(_node, _triangle) {}
    
    Float operator() ( const BVHBox& bbox_left, const int triangle_left, const BVHBox& bbox_right, const int triangle_right ) const
    {
        Float nodeSA= BVHBox(bbox_left, bbox_right).area();
        
        return node_cost * 2
            + bbox_left.area() / nodeSA * triangle_left * triangle_cost
            + bbox_right.area() / nodeSA * triangle_right * triangle_cost;
    }
    
    Float operator() ( const int triangles ) const
    {
        return triangles * triangle_cost;
    }
};


struct compareBucket
{
    BVHBox bounds;
    int split, count, axis;
    Float scale;
    
    compareBucket( const int _split, const int _n, const int _axis, const BVHBox& _bbox )
        : bounds(_bbox), split(_split), count(_n), axis(_axis),  scale( _n / (_bbox.pmax(_axis) - _bbox.pmin(_axis)) ) {}
    
    bool operator() ( const BVHPrimitive& primitive ) const
    {
        int b= (primitive.center(axis) - bounds.pmin(axis)) * scale;
        if(b < 0) b= 0;
        if(b > count -1) b= count -1;
        return (b <= split);
    }
};


struct BVHBucket
{
    BVHBox bounds;
    int count;
    
    BVHBucket( ) : bounds(), count(0) {}
};


/*! standard binning algorithm using a cost (SAH) function.. 
    wald 2007
    "On fast Construction of SAH-based Bounding Volume Hierarchies"
    http://www.sci.utah.edu/~wald/Publications/2007/ParallelBVHBuild/fastbuild.pdf    
*/
class BVHBuilder
{
public:
    BVHBuilder( const BVHCost& cost, const unsigned int max_primitives= 1 );
    ~BVHBuilder( );
    
    void build( BVH& _bvh ); //!< attention : les triangles doivent etre affectes au bvh avant !!
    
protected:
    std::vector<BVHNode> nodes;
    std::vector<BVHPrimitive> primitives;
    std::vector<int> node_axis;

    const BVHCost& eval;
    int max_primitives;
    
    //! internal build procedure.
    int build_node( const int begin, const int end );
};

#endif

