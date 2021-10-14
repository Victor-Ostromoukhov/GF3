
#include <cmath>
#include <algorithm>

#include "bvh_builder.h"
#include "stbvh_builder.h"


// \todo version simplifiee : construit un englobant pour un intervalle exact t0 et t1... 
// interpolation supplementaire si l'intervalle r0 r1 n'est pas exactement t0 t1, mais t0 t1+1, pour t= 0.. 1
STBVHPrimitive STBVHBuilder::build_primitive( const STBVHMesh& object, const int index, const int t0, const int t1 )
{
    assert(t0 < int(object.vertices.size()));
    assert(t1 < int(object.vertices.size()));
    assert(t0 <= t1);
    assert(object.vertices[t0]->size() == object.vertices[t1]->size());

    int v0= object.indices[index];
    int v1= object.indices[index +1];
    int v2= object.indices[index +2];
    
    BVHBox b0;
    b0.insert(Point((*object.vertices[t0])[v0]));
    b0.insert(Point((*object.vertices[t0])[v1]));
    b0.insert(Point((*object.vertices[t0])[v2]));
    
    BVHBox b1;
    b1.insert(Point((*object.vertices[t1])[v0]));
    b1.insert(Point((*object.vertices[t1])[v1]));
    b1.insert(Point((*object.vertices[t1])[v2]));
    
    for(int i= 0; i < 3; i++)
    {
        Float mtmin= 0;
        Float mtmax= 0;
        
        for(int k= t0; k <= t1; k++)
        {
            Float t= Float(k - t0) / Float(t1 - t0);
            if(t < 0) t= 0;
            if(t > 1) t= 1;
            Float tmin= b0.pmin(i) * (1 - t) + b1.pmin(i) * t;
            Float tmax= b0.pmax(i) * (1 - t) + b1.pmax(i) * t;
            
            // test : no espilon            
            if((*object.vertices[k])[v0](i) < tmin) { mtmin= std::max(mtmin, tmin - (*object.vertices[k])[v0](i) ); }
            if((*object.vertices[k])[v1](i) < tmin) { mtmin= std::max(mtmin, tmin - (*object.vertices[k])[v1](i) ); }
            if((*object.vertices[k])[v2](i) < tmin) { mtmin= std::max(mtmin, tmin - (*object.vertices[k])[v2](i) ); }
            
            if((*object.vertices[k])[v0](i) > tmax) { mtmax= std::max(mtmax, (*object.vertices[k])[v0](i) - tmax); }
            if((*object.vertices[k])[v1](i) > tmax) { mtmax= std::max(mtmax, (*object.vertices[k])[v1](i) - tmax); }
            if((*object.vertices[k])[v2](i) > tmax) { mtmax= std::max(mtmax, (*object.vertices[k])[v2](i) - tmax); }
        }
        
        b0.pmin(i)= b0.pmin(i) - mtmin;
        b1.pmin(i)= b1.pmin(i) - mtmin;
        b0.pmax(i)= b0.pmax(i) + mtmax;
        b1.pmax(i)= b1.pmax(i) + mtmax;
    }
    

    #if 0    
    // sanity check
    {
        for(int i= 0; i < 3; i++)
        {
            for(int k= t0; k <= t1; k++)
            {
                Float t= Float(k - t0) / Float(t1 - t0);
                if(t < 0) t= 0;
                if(t > 1) t= 1;
                Float tmin= b0.pmin(i) * (1 - t) + b1.pmin(i) * t;
                Float tmax= b0.pmax(i) * (1 - t) + b1.pmax(i) * t;
                
                if((*object.vertices[k])[v0](i) < tmin) printf("index %d: min %f %f = %.10f\n", index, (*object.vertices[k])[v0](i), tmin, tmin - (*object.vertices[k])[v0](i));
                if((*object.vertices[k])[v1](i) < tmin) printf("index %d: min %f %f = %.10f\n", index, (*object.vertices[k])[v1](i), tmin, tmin - (*object.vertices[k])[v1](i));
                if((*object.vertices[k])[v2](i) < tmin) printf("index %d: min %f %f = %.10f\n", index, (*object.vertices[k])[v2](i), tmin, tmin - (*object.vertices[k])[v2](i));
                
                if((*object.vertices[k])[v0](i) > tmax) printf("index %d: max %f %f = %.10f\n", index, (*object.vertices[k])[v0](i), tmax, (*object.vertices[k])[v0](i) - tmax);
                if((*object.vertices[k])[v1](i) > tmax) printf("index %d: max %f %f = %.10f\n", index, (*object.vertices[k])[v1](i), tmax, (*object.vertices[k])[v1](i) - tmax);
                if((*object.vertices[k])[v2](i) > tmax) printf("index %d: max %f %f = %.10f\n", index, (*object.vertices[k])[v2](i), tmax, (*object.vertices[k])[v2](i) - tmax);
            }
        }
    }
    #endif
    
    return STBVHPrimitive(b0, b1, object.index, index, object.slots, t0, t1);
}


void STBVHBuilder::build_bounds( const STBVH& bvh, std::vector<STBVHPrimitive>& primitives, const Float r0, const Float r1, BVHBox& b0, BVHBox& b1 )
{
    b0.clear();
    b1.clear();
    
    for(int i= 0; i < int(primitives.size()); i++)
    {
        //~ if(primitives[i].motion(r0, r1) == 0)
            //~ continue;
        
        b0.insert(primitives[i].bounds0);
        b1.insert(primitives[i].bounds1);
    }
    
    for(int axis= 0; axis < 3; axis++)
    {
        Float mtmin= 0;
        Float mtmax= 0;
        
        for(int i= 0; i < int(primitives.size()); i++)
        {
            //~ if(primitives[i].motion(r0, r1) == 0)
                //~ continue;
            
            assert(primitives[i].object < int(bvh.objects.size()));
            const std::vector<unsigned int>& indices= bvh.objects[primitives[i].object].indices;
            const std::vector< const std::vector<vec3> * >& vertices= bvh.objects[primitives[i].object].vertices;            
            
            int index= primitives[i].index;
            assert(index + 2 < int(indices.size()));
            int v0= indices[index];
            int v1= indices[index +1];
            int v2= indices[index +2];
            
            int t0= primitives[i].slot0(r0);
            int t1= primitives[i].slot1(r1);
            assert(t0 < int(vertices.size()));
            assert(t1 < int(vertices.size()));
            
            for(int k= t0; k <= t1; k++)
            {
                Float t= Float(k - t0) / Float(t1 - t0);
                if(t < 0) t= 0;
                if(t > 1) t= 1;
                Float tmin= b0.pmin(axis) * (1 - t) + b1.pmin(axis) * t;
                Float tmax= b0.pmax(axis) * (1 - t) + b1.pmax(axis) * t;
                
                // test : no espilon
                if((*vertices[k])[v0](axis) < tmin) { mtmin= std::max( mtmin, tmin - (*vertices[k])[v0](axis) ); }
                if((*vertices[k])[v1](axis) < tmin) { mtmin= std::max( mtmin, tmin - (*vertices[k])[v1](axis) ); }
                if((*vertices[k])[v2](axis) < tmin) { mtmin= std::max( mtmin, tmin - (*vertices[k])[v2](axis) ); }
                
                if((*vertices[k])[v0](axis) > tmax) { mtmax= std::max( mtmax, (*vertices[k])[v0](axis) - tmax ); }
                if((*vertices[k])[v1](axis) > tmax) { mtmax= std::max( mtmax, (*vertices[k])[v1](axis) - tmax ); }
                if((*vertices[k])[v2](axis) > tmax) { mtmax= std::max( mtmax, (*vertices[k])[v2](axis) - tmax ); }
            }
        }
        
        b0.pmin(axis)= b0.pmin(axis) - mtmin;
        b1.pmin(axis)= b1.pmin(axis) - mtmin;
        b0.pmax(axis)= b0.pmax(axis) + mtmax;
        b1.pmax(axis)= b1.pmax(axis) + mtmax;
    }

    #if 0    
    // sanity check
    {
        for(int i= 0; i < 3; i++)
        {
            for(int k= t0; k <= t1; k++)
            {
                Float t= Float(k - t0) / Float(t1 - t0);
                Float tmin= b0.pmin(i) * (1 - t) + b1.pmin(i) * t;
                Float tmax= b0.pmax(i) * (1 - t) + b1.pmax(i) * t;
                    
                for(int j= 0; j < int(primitives.size()); j++)
                {
                    int index= primitives[j].index;
                    int v0= bvh.indices[index];
                    int v1= bvh.indices[index +1];
                    int v2= bvh.indices[index +2];
                    
                    if((*bvh.vertices[k])[v0](i) < tmin) printf("bounds %d: min %f %f = %.10f\n", index, (*bvh.vertices[k])[v0](i), tmin, tmin - (*bvh.vertices[k])[v0](i));
                    if((*bvh.vertices[k])[v1](i) < tmin) printf("bounds %d: min %f %f = %.10f\n", index, (*bvh.vertices[k])[v1](i), tmin, tmin - (*bvh.vertices[k])[v1](i));
                    if((*bvh.vertices[k])[v2](i) < tmin) printf("bounds %d: min %f %f = %.10f\n", index, (*bvh.vertices[k])[v2](i), tmin, tmin - (*bvh.vertices[k])[v2](i));
                    
                    if((*bvh.vertices[k])[v0](i) > tmax) printf("bounds %d: max %f %f = %.10f\n", index, (*bvh.vertices[k])[v0](i), tmax, (*bvh.vertices[k])[v0](i) - tmax);
                    if((*bvh.vertices[k])[v1](i) > tmax) printf("bounds %d: max %f %f = %.10f\n", index, (*bvh.vertices[k])[v1](i), tmax, (*bvh.vertices[k])[v1](i) - tmax);
                    if((*bvh.vertices[k])[v2](i) > tmax) printf("bounds %d: max %f %f = %.10f\n", index, (*bvh.vertices[k])[v2](i), tmax, (*bvh.vertices[k])[v2](i) - tmax);
                }
            }
        }
    }
    #endif
}


void STBVHBuilder::build( STBVH& bvh )
{
    bvh.nodes.clear();
    bvh.references.clear();
    bvh.bounds.clear();
    
    // static tree ?
    bvh.notime= true;
    for(int i= 0; i < int(bvh.objects.size()); i++)
        if(bvh.objects[i].buffers() > 1)
            bvh.notime= false;
    
    // static objects : force 1 segment
    for(int i= 0; i < int(bvh.objects.size()); i++)
        if(bvh.objects[i].buffers() == 1)
        {
            bvh.objects[i].buffers(2);
            bvh.objects[i].vertices[1]= bvh.objects[i].vertices[0];
        }
        
    int n= 0;
    for(int i= 0; i < int(bvh.objects.size()); i++)
        n= n + int(bvh.objects[i].indices.size()) / 3;
    
    std::vector<STBVHPrimitive> primitives;
    primitives.reserve(n);
    references.reserve(n * 2);
    
    for(int k= 0; k < int(bvh.objects.size()); k++)
        for(int i= 0; i < int(bvh.objects[k].indices.size()); i+= 3)
        {
            assert(bvh.objects[k].vertices.size() > 1);
            primitives.push_back( build_primitive(bvh.objects[k], i, /* t0 */ 0, /* t1 */ int(bvh.objects[k].vertices.size()) -1) );
        }
        
    bvh.root= build_node(bvh, primitives, /* r0 */ 0, /* r1 */ 1);
    std::swap(bvh.nodes, nodes);
    std::swap(bvh.references, references);
    
    int nodes= 0;
    int leafs= 0;
    for(int i= 0; i < int(bvh.nodes.size()); i++)
        if(bvh.nodes[i].internal())
            nodes++;
        else
            leafs++;
    
    printf("\nlinear motion BVH %d nodes: %d internal, %d leafs, root %d\n", int(bvh.nodes.size()), nodes, leafs, bvh.root);
    printf("  nodes %d %dMB\n", int(bvh.nodes.size()), int(bvh.nodes.size() * sizeof(STBVHNode) / 1024 / 1024));
    printf("    internals %d %dMB\n", nodes, int(nodes * sizeof(STBVHNode) / 1024 / 1024));
    printf("    leafs %d %dMB\n", leafs, int(leafs * sizeof(STBVHNode) / 1024 / 1024));
    printf("  references %dMB\n", int(bvh.references.size() * sizeof(int) / 1024 / 1024));
    
    int buffers= 0;
    size_t buffer_size= 0;
    size_t indice_size= 0;
    for(int k= 0; k < int(bvh.objects.size()); k++)
    {
        indice_size= indice_size + bvh.objects[k].indices.size() * sizeof(unsigned int);
        buffers= buffers + int(bvh.objects[k].vertices.size());
        for(int i= 0; i < int(bvh.objects[k].vertices.size()); i++)
            buffer_size= buffer_size + bvh.objects[k].vertices[i]->size() * sizeof(vec3);
    }
    
    printf("  buffers %d %dMB\n", buffers, int(buffer_size / 1024 / 1024));
    printf("  indices %dMB\n", int(indice_size / 1024 / 1024));
    
    if(bvh.notime) 
        printf("  all objects are static...\n");
    printf("\n");
}


struct compareSTBucket
{
    BVHBox bounds;
    int split, count, axis;
    Float scale;
    
    compareSTBucket( const int _split, const int _n, const int _axis, const BVHBox& _bbox )
        : bounds(_bbox), split(_split), count(_n), axis(_axis),  scale( _n / (_bbox.pmax(_axis) - _bbox.pmin(_axis)) ) {}
    
    bool operator() ( const STBVHPrimitive& primitive ) const
    {
        int b= (primitive.center05(axis) - bounds.pmin(axis)) * scale;
        if(b < 0) b= 0;
        if(b > count -1) b= count -1;
        return (b <= split);
    }
};


int STBVHBuilder::build_node( const STBVH& bvh, std::vector<STBVHPrimitive>& primitives, const Float r0, const Float r1 )
{
    if(int(primitives.size()) <= 2)
    {
        BVHBox node_b0, node_b1;
        build_bounds(bvh, primitives, r0, r1, node_b0, node_b1);
        
        int first= int(references.size());
        for(int i= 0; i < int(primitives.size()); i++)
            references.emplace_back(primitives[i].object, primitives[i].index);
        
        //~ printf("force leaf %d: time %f %f, %d triangles\n", int(nodes.size()), r0, r1, int(references.size()) - first);
        
        // construction feuille
        STBVHNode node;
        node.make_leaf(node_b0, r0, node_b1, r1, first, int(references.size()));
        nodes.push_back(node);
        return int(nodes.size()) -1;
    }

    // leaf cost
    BVHBox bounds;
    BVHBox centroids;
    Float leaf_cost= 0;
    for(int i= 0; i < int(primitives.size()); i++)
    {
        //~ if(primitives[i].motion(r0, r1) == 0)
            //~ continue;

        centroids.insert(primitives[i].center05);
        bounds.insert(primitives[i].bounds05);
        leaf_cost= leaf_cost + primitives[i].motion(r0, r1);
    }
    
    // object cost
    int minAxis= -1;
    int minCostSplit= -1;
    Float minCost= FLOAT_MAX;
    const int nBuckets= 16;
    for(int axis= 0; axis < 3; axis++)
    {
        BVHBucket buckets[nBuckets];
        
        Float scale= Float(nBuckets) / (centroids.pmax(axis) - centroids.pmin(axis));
        //~ if(std::isinf(scale))
            //~ continue;
        
        for(int i= 0; i < int(primitives.size()); i++)
        {
            int b= (primitives[i].center05(axis) - centroids.pmin(axis)) * scale;
            if(b < 0) b= 0;
            if(b > nBuckets -1) b= nBuckets -1;
            
            buckets[b].bounds.insert(primitives[i].bounds05);
            buckets[b].count+= primitives[i].motion(r0, r1);
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
            
            // \todo merge left/right sweep !!
            
            Float cost= 2 + b0.area() / bounds.area() * count0 + b1.area() / bounds.area() * count1;
            // object split T(Y) / T(X) == 1 
            if(cost < minCost)
            {
                minAxis= axis;
                minCost= cost;
                minCostSplit= i;
            }
        }
    }
    //~ assert(minAxis != -1);
    //~ assert(minCostSplit >= 0);
    Float object_cost= minCost;
    
    // temporal cost
    Float temporal_cost= FLOAT_MAX;
    
    // round split time
    int motion= 0;
    int motion_max= -1;
    for(int i= 0; i < int(primitives.size()); i++)
    {
        int m= primitives[i].motion(r0, r1);
        if(m > motion)
        {
            motion= m;
            motion_max= i;
        }
    }
    assert(motion_max != -1);
    Float r05= primitives[motion_max].round0((r0 + r1) / 2);
    
    Float r00= 1;   // time range 0 / left
    Float r01= 0;
    Float r10= 1;   // time range 1 / right
    Float r11= 0;
    
    if(object_cost >= leaf_cost / 2 && motion > 1 )
    {
        // evalue la coupe temporelle si necessaire, reduit la taille de l'arbre...
        BVHBox b0, b1;
        int count0= 0, count1= 0;
        
        for(int i= 0; i < int(primitives.size()); i++)
        {
            // time range T(X0) 
            r00= std::min(r00, primitives[i].round0(r0));
            r01= std::max(r01, primitives[i].round1(r05));
            
            // time range T(X1)
            r10= std::min(r10, primitives[i].round0(r05));
            r11= std::max(r11, primitives[i].round1(r1));
            
            if(primitives[i].motion(r0, r1) > 1)
            {
                STBVHPrimitive p0= build_primitive(bvh.objects[primitives[i].object], primitives[i].index, primitives[i].slot0(r0), primitives[i].slot1(r05));
                b0.insert(p0.bounds05);
                count0+= primitives[i].motion(r0, r05);
                
                STBVHPrimitive p1= build_primitive(bvh.objects[primitives[i].object], primitives[i].index, primitives[i].slot0(r05), primitives[i].slot1(r1));
                b1.insert(p1.bounds05);
                count1+= primitives[i].motion(r05, r1);
            }
            else // if(primitives[i].motion(r0, r1) == 1)
            {
                // objet statique, dupliquer la primitive...
                b0.insert(primitives[i].bounds0);
                b0.insert(primitives[i].bounds1);
                count0+= primitives[i].motion(r0, r1);
                
                b1.insert(primitives[i].bounds0);
                b1.insert(primitives[i].bounds1);
                count1+= primitives[i].motion(r0, r1);
            }
        }
        
        // temporal split T(Y) / T(X) 
        temporal_cost= 2 + b0.area() / bounds.area() * count0 * (r01 - r00) / (r1 - r0) + b1.area() / bounds.area() * count1 * (r11 - r10) / (r1 - r0) ;
        temporal_cost= Float(1.25) * temporal_cost;
    }

    //~ if(temporal_cost < FLOAT_MAX)
        //~ printf("node %d %d: time %f %f: leaf %f, object %f, temporal %d %f\n", int(nodes.size()), int(primitives.size()), r0, r1, leaf_cost, object_cost, motion, temporal_cost);
    
    //
    BVHBox node_b0, node_b1;
    build_bounds(bvh, primitives, r0, r1, node_b0, node_b1);
    
    if(leaf_cost < object_cost && leaf_cost < temporal_cost)
    {
        int first= int(references.size());
        for(int i= 0; i < int(primitives.size()); i++)
            references.emplace_back(primitives[i].object, primitives[i].index);
        
        //~ printf("  leaf %d: time %f %f, %d triangles\n", int(nodes.size()), r0, r1, int(references.size()) - first);
        
        // construction feuille
        STBVHNode node;
        node.make_leaf(node_b0, r0, node_b1, r1, first, int(references.size()));
        nodes.push_back(node);
        return int(nodes.size()) -1;
    }
    
    int left= -1;
    int right= -1;
    if(object_cost == FLOAT_MAX || object_cost < temporal_cost)
    {
        //~ assert(minAxis != -1);
        //~ assert(minCostSplit >= 0);
        if(minAxis == -1)
        {
            minAxis= 0;
            minCostSplit= nBuckets / 2;
            printf("[stbvh] min axis == -1, %d primitives\n", int(primitives.size()));
        }
        
        // object partition 
        STBVHPrimitive *p= std::partition( primitives.data(), primitives.data() + primitives.size(), compareSTBucket(minCostSplit, nBuckets, minAxis, centroids) );
        int m= int(std::distance(primitives.data(), p));
        
        if(m == 0 || m == int(primitives.size()))
            // si la partion est degeneree, force un median split
            m= int(primitives.size()) / 2;
        assert(m != 0);
        assert(m != int(primitives.size()));
    
        //~ printf("  object %d: time %f %f, %d / %d triangles\n", int(nodes.size()), r0, r1, int(primitives.size()), m);
        
        std::vector<STBVHPrimitive> tmp;
        tmp.reserve(primitives.size() - m);
        for(int i= m; i < int(primitives.size()); i++)
            tmp.push_back(primitives[i]);
        
        primitives.resize(m);
        
        left= build_node(bvh, primitives, r0, r1);
        right= build_node(bvh, tmp, r0, r1);
        
        // \todo merge bounds
    }
    else 
    {
        assert(temporal_cost < FLOAT_MAX);
        
        // temporal partition
        std::vector<STBVHPrimitive> tmp;
        tmp.reserve(primitives.size());
        for(int i= 0; i < int(primitives.size()); i++)
            if(primitives[i].motion(r0, r1) > 1)
                tmp.push_back( build_primitive(bvh.objects[primitives[i].object], primitives[i].index, primitives[i].slot0(r05), primitives[i].slot1(r1)) );
            else
                tmp.push_back( primitives[i] );
            
        for(int i= 0; i < int(primitives.size()); i++)
            if(primitives[i].motion(r0, r1) > 1)
                primitives[i]= build_primitive(bvh.objects[primitives[i].object], primitives[i].index, primitives[i].slot0(r0), primitives[i].slot1(r05));
            //~ else
                //~ primitives[i]= primitives[i];
            
        //~ printf("  temporal %d: time %f %f / %f, %d triangles\n", int(nodes.size()), r0, r1, r05, int(primitives.size()));
        //~ printf("    time %f %f / %f %f\n", r00, r01, r10, r11);
        
        left= build_node(bvh, primitives, r0, r05);
        right= build_node(bvh, tmp, r05, r1);
        
        // \todo merge bounds
    }
    
    // \todo memory !!
    
    STBVHNode node;
    node.make_node(node_b0, std::min(nodes[left].r0, nodes[right].r0), node_b1, std::max(nodes[left].r1, nodes[right].r1), left, right);
    nodes.push_back(node);
    return int(nodes.size()) -1;
}

