
#include <cassert>
#include <algorithm>

#include "bvh.h"


// predicat pour repartir les triangles avant / apres la position de la coupe choisie pour separer les fils d'un noeud
struct before
{
    int axis;
    Float cut;
    
    before( const int _axis, const Float _cut ) : axis(_axis), cut(_cut) {}
    
    bool operator( ) ( const BVHPrimitive& primitive ) const
    {
        return primitive.center(axis) <= cut;
    }
};


// construction recursive d'un noeud
int build_node( std::vector<BVHNode>& nodes, std::vector<BVHPrimitive>& primitives, const int begin, const int end )
{
    // feuille
    if(begin +1 >= end)
    {
        BVHNode node; 
        node.make_leaf( primitives[begin].bounds, begin, end );
        nodes.push_back(node);
        return (int) nodes.size() - 1;
    }
    
    // englobant des centres des objets
    BVHBox cbox;
    for(int i= begin; i < end; i++)
        cbox.insert(primitives[i].center);
    
    // axe et coupe
    int axis= 0;
    Vector d= cbox.pmax - cbox.pmin;
    if(d.x > d.y && d.x > d.z) axis= 0;
    else if(d.y > d.x && d.y > d.z) axis= 1;
    else axis= 2;
    
    Float cut= cbox.pmin(axis) + d(axis) / 2;
    
    // partition
    BVHPrimitive *pmid= std::partition(primitives.data() + begin, primitives.data() + end, before(axis, cut));
    int mid= std::distance(primitives.data(), pmid);
    
    // cas degenere, reparti arbitrairement les objets en 2
    if(mid == begin || mid == end)
        mid= (begin + end) / 2;
    assert(mid != begin);
    assert(mid != end);
    
    // construction recursive des fils
    int left= build_node(nodes, primitives, begin, mid);
    int right= build_node(nodes, primitives, mid, end);

    // construction du noeud
    BVHNode node;
    node.make_node( BVHBox(nodes[left].bounds, nodes[right].bounds), left, right );
    nodes.push_back(node);
    return (int) nodes.size() -1;
}


BVH::BVH( ) {}

void BVH::build( )
{
    assert(triangles.size());
    
    std::vector<BVHPrimitive> primitives;
    primitives.reserve(triangles.size());
    for(unsigned int i= 0; i < triangles.size(); i++)
        primitives.push_back( BVHPrimitive(triangles[i].bounds(), i) );
    
    root= build_node(nodes, primitives, 0, (int) primitives.size());
    bounds= nodes[root].bounds;
    
    std::vector<BVHTriangle> remap;
    remap.reserve(primitives.size());
    for(unsigned int i= 0; i < primitives.size(); i++)
        remap.push_back( BVHTriangle(triangles[primitives[i].object_id], primitives[i].object_id) );
    triangles.swap(remap);
}


void BVH::intersect( const int index, const Ray& ray, const Vector& invd, Hit& hit ) const
{
    if(index == -1)
        return;
    
    const BVHNode& node= nodes[index];
    if(node.internal())
    {
        // parcours ordonne
        Float left_min, left_max;
        bool left= nodes[node.left].bounds.intersect(ray, invd, hit.t, left_min, left_max);

        Float right_min, right_max;
        bool right= nodes[node.right].bounds.intersect(ray, invd, hit.t, right_min, right_max);
        
        // parcours ordonne
        if(left && right)
        {
            Float left_mid= (left_min + left_max) / 2;
            Float right_mid= (right_min + right_max) / 2;
            if(left_mid < right_mid)
            {
                intersect(node.left, ray, invd, hit);
                if(hit.t >= right_min)  // ne teste l'autre fils que si necessaire...
                    intersect(node.right, ray, invd, hit);
            }
            else
            {
                intersect(node.right, ray, invd, hit);
                if(hit.t >= left_min)
                    intersect(node.left, ray, invd, hit);
            }
        }
        else if(left)
            intersect(node.left, ray, invd, hit);
        else if(right)
            intersect(node.right, ray, invd, hit);
    }
    
    else // if(node.leaf())
    {
        int begin= node.leaf_begin();
        int end= node.leaf_end();
        for(int i= begin; i < end; i++)
        {
            Float t, u, v;
            if(triangles[i].intersect(ray, hit.t, t, u, v))
            {
                hit.t= t;   // !! raccourcir le rayon !! evite de visiter trop de noeuds apres avoir trouve une intersection
                hit.u= u;
                hit.v= v;
                hit.triangle_id=triangles[i].id;   // indice du triangle dans l'ensemble initial
            }
        }
    }
}


bool BVH::intersect( const Ray& ray, Hit& hit ) const
{
    Vector invd(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);

    hit.t= ray.tmax;
    hit.triangle_id= -1;
    
#if 0
    // version recursive
    intersect(root, ray, invd, hit);
#else

    // version "manuelle"
    int stack[64];
    int top= 0;
    
    int index= root;
    while(true)
    {
        const BVHNode& node= nodes[index];
        if(node.internal())
        {
            // parcours ordonne
            Float left_min, left_max;
            bool left= nodes[node.left].bounds.intersect(ray, invd, hit.t, left_min, left_max);
            
            Float right_min, right_max;
            bool right= nodes[node.right].bounds.intersect(ray, invd, hit.t, right_min, right_max);
            
            // parcours ordonne
            if(left && right)
            {
                int push;
                Float left_mid= (left_min + left_max) / 2;
                Float right_mid= (right_min + right_max) / 2;
                if(left_mid < right_mid)
                {
                    index= node.left;
                    push= node.right;
                }
                else
                {
                    index= node.right;
                    push= node.left;
                }
                
                // push
                stack[top++]= push;
            }
            else if(left)
                index= node.left;
            else if(right)
                index= node.right;
            else
            {
                // pop
                if(top == 0) break;
                index= stack[--top];
            }
        }
        
        else // if(node.leaf())
        {
            Float t, u, v;
            for(int i= node.leaf_begin(); i < node.leaf_end(); i++)
            {
                if(triangles[i].intersect(ray, hit.t, t, u, v))
                {
                    hit.t= t;   // !! raccourcir le rayon !! evite de visiter trop de noeuds apres avoir trouve une intersection
                    hit.u= u;
                    hit.v= v;
                    hit.triangle_id=triangles[i].id;   // indice du triangle dans l'ensemble initial
                }
            }
            
            // pop
            if(top == 0) break;
            index= stack[--top];
        }
    }
#endif
    
    return (hit.triangle_id != -1);
}


bool BVH::visible( const int index, const Ray& ray, const Vector& invd ) const
{
    const BVHNode& node= nodes[index];
    if(node.internal())
    {
        Float tmin, tmax;
        if(node.bounds.intersect(ray, invd, ray.tmax, tmin, tmax)) 
        {
            if(!visible(node.left, ray, invd))
                return false;
            if(!visible(node.right, ray, invd))
                return false;
        }
    }

    else // if(node.leaf())
    {
        Float t, u, v;
        for(int i= node.leaf_begin(); i < node.leaf_end(); i++)
            if(triangles[i].intersect(ray, ray.tmax, t, u, v))
                return false;
    }
    
    return true;
}

bool BVH::visible( const Point& p, const Point& q ) const
{
    assert(root != -1);
    
    Ray ray(p, q);
    Vector invd(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    
#if 0
    // version recursive
    return visible(root, ray, invd);
#else

    // version "manuelle"
    int stack[64];
    int top= 0;
    
    int index= root;
    while(true)
    {
        const BVHNode& node= nodes[index];
        if(node.internal())
        {
            Float left_min, left_max;
            bool left= nodes[node.left].bounds.intersect(ray, invd, ray.tmax, left_min, left_max);
            
            Float right_min, right_max;
            bool right= nodes[node.right].bounds.intersect(ray, invd, ray.tmax, right_min, right_max);
            
            if(left && right)
            {
            #if 0
                // parcours ordonne
                int push;
                Float left_mid= (left_min + left_max) / 2;
                Float right_mid= (right_min + right_max) / 2;
                if(left_mid < right_mid)
                {
                    index= node.left;
                    push= node.right;
                }
                else
                {
                    index= node.right;
                    push= node.left;
                }
                
                // push
                stack[top++]= push;
            #else
            
                // parcours prefixe, plus rapide ??
                stack[top++]= node.right;
                index= node.left;
            #endif
            }
            else if(left)
                index= node.left;
            else if(right)
                index= node.right;
            else
            {
                // pop
                if(top == 0) 
                    return true;
                index= stack[--top];
            }
        }
        
        else // if(node.leaf())
        {
            Float t, u, v;
            for(int i= node.leaf_begin(); i < node.leaf_end(); i++)
                if(triangles[i].intersect(ray, ray.tmax, t, u, v))
                    return false;
            
            // pop
            if(top == 0) 
                return true;
            index= stack[--top];
        }
    }
#endif
}

