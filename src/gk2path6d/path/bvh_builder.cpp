
#include <algorithm>

#include "bvh_builder.h"


BVHBuilder::BVHBuilder( const BVHCost &_eval, const unsigned int _max_primitives ) 
    : nodes(), primitives(), eval(_eval), max_primitives(_max_primitives) {}
BVHBuilder::~BVHBuilder( ) {}

    
int BVHBuilder::build_node( const int begin, const int end )
{
    // force la creation d'une feuille si 1 primitive, laisse le sah decider de la taille des feuilles
    if(begin + 1 >= end)
    {
        // calcule la boite englobante du sous ensemble de triangles
        BVHBox bbox;
        for(int i= begin; i < end; i++)
            bbox.insert(primitives[i].bounds);
        
        BVHNode node;
        node.make_leaf(bbox, begin, end);
        nodes.push_back(node);
        node_axis.push_back(-1);
        return (int) nodes.size() -1;
    }
    
    // calcule la boite englobante des centres du sous ensemble de triangles
    BVHBox centroids;
    for(int i= begin; i < end; i++)
        centroids.insert(primitives[i].center);
    
    // cf pbrt 2
    int minAxis= -1;
    int minCostSplit= -1;
    Float minCost= FLOAT_MAX;
    const int nBuckets= 16;
    for(int axis= 0; axis < 3; axis++)
    {
        BVHBucket buckets[nBuckets];
        
        Float scale= Float(nBuckets) / (centroids.pmax(axis) - centroids.pmin(axis));
        for(int i= begin; i < end; i++)
        {
            int b= (primitives[i].center(axis) - centroids.pmin(axis)) * scale;
            if(b < 0) b= 0;
            if(b > nBuckets -1) b= nBuckets -1;
            
            buckets[b].bounds.insert(primitives[i].bounds);
            buckets[b].count++;
        }
        
        for(int i= 0; i < nBuckets -1; i++ )
        {
            BVHBox b0, b1;
            int count0= 0, count1= 0;
            
            for(int j= 0; j <= i; j++)     // compte les triangles avant la coupe...
            {
                b0.insert(buckets[j].bounds);
                count0= count0 + buckets[j].count;
            }
            
            for(int j= i +1; j < nBuckets; j++)   // ... et les triangles restant, apres la coupe
            {
                b1.insert(buckets[j].bounds);
                count1= count1 + buckets[j].count;
            }
            
            // Find bucket to split at that minimizes SAH metric
            Float cost= eval(b0, count0, b1, count1);
            if(cost < minCost)
            {
                minAxis= axis;
                minCost= cost;
                minCostSplit= i;
            }
        }
    }
    assert(minAxis != -1);
    assert(minCostSplit >= 0);
    
    // force la construction des feuilles, si la coupe n'est pas interressante
    if(end - begin < 16 && minCost > eval(end - begin))
    {
        BVHBox bbox;
        for(int i= begin; i < end; i++)
            bbox.insert(primitives[i].bounds);
        
        BVHNode node;
        node.make_leaf(bbox, begin, end);
        
        assert(minAxis != -1);
        assert(node_axis.size() == nodes.size());
        nodes.push_back(node);
        node_axis.push_back(minAxis);
        return (int) nodes.size() -1;
    }

    BVHPrimitive *p= std::partition( primitives.data() + begin, primitives.data() + end, compareBucket(minCostSplit, nBuckets, minAxis, centroids) );
    int m= (int) std::distance(primitives.data(), p);
    
    if(m == begin || m == end)
        // si la partion est degeneree, force un median split
        m= ( begin + end ) / 2;
    assert(m != begin);
    assert(m != end);
    
    // construction recursive
    int left= build_node(begin, m);
    int right= build_node(m, end);
    
    // construire le noeud interne
    BVHNode node;
    node.make_node(BVHBox(nodes[left].bounds, nodes[right].bounds), left, right);
    
    // stocke le noeud
    assert(minAxis != -1);
    assert(node_axis.size() == nodes.size());
    nodes.push_back(node);
    node_axis.push_back(minAxis);
    
    return (int) nodes.size() -1;
}


void BVHBuilder::build( BVH& bvh )
{
    assert(bvh.triangles.size());
    
    bvh.bounds.clear();
    bvh.nodes.clear();
    bvh.node_axis.clear();
    bvh.root= -1;
    
    primitives.reserve(bvh.triangles.size());
    for(unsigned int i= 0; i < bvh.triangles.size(); i++)
    {
        BVHBox box= bvh.triangles[i].bounds();
        primitives.push_back( BVHPrimitive(box, i) );
    }
    
    // build
    nodes.reserve(primitives.size() * 2 / max_primitives);
    bvh.root= build_node(0, (int) primitives.size());
    bvh.bounds= nodes[bvh.root].bounds;
    
    assert(node_axis.size() == nodes.size());
    bvh.nodes.swap(nodes);
    bvh.node_axis.swap(node_axis);
    
    printf("BVHBuilder: %d nodes, root %d\n", (int) bvh.nodes.size(), bvh.root);
    
    // reorder triangles
    std::vector<BVHTriangle> remap;
    remap.reserve(primitives.size());
    for(unsigned int i= 0; i < primitives.size(); i++)
        remap.push_back( bvh.triangles[primitives[i].object_id] );
    
    bvh.triangles.swap(remap);
}

