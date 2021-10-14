
#include <cmath>

#include "stbvh.h"


void STBVH::intersect( const int index, const Ray& ray, const Vector& invd, Hit& hit ) const
{
    assert(index != -1);
    
    const STBVHNode& node= nodes[index];
    if(ray.time < node.r0 || ray.time > node.r1) return;
    
    if(node.internal())
    {
        // interpole les bbox des fils
        Float tleft= (ray.time - nodes[node.left].r0) / (nodes[node.left].r1 - nodes[node.left].r0);
        Float left_min, left_max;
        bool left= false;
        if(tleft >= 0 && tleft <= 1)
        {
            BVHBox left_bounds= BVHBox(nodes[node.left].bounds0.pmin * (1 - tleft) + nodes[node.left].bounds1.pmin * tleft, nodes[node.left].bounds0.pmax * (1 - tleft) + nodes[node.left].bounds1.pmax * tleft);
            left= left_bounds.intersect(ray, invd, hit.t, left_min, left_max);
        }
        
        Float tright= (ray.time - nodes[node.right].r0) / (nodes[node.right].r1 - nodes[node.right].r0);
        Float right_min, right_max;
        bool right= false;
        if(tright >= 0 && tright <= 1)
        {
            BVHBox right_bounds= BVHBox(nodes[node.right].bounds0.pmin * (1 - tright) + nodes[node.right].bounds1.pmin * tright, nodes[node.right].bounds0.pmax * (1 - tright) + nodes[node.right].bounds1.pmax * tright);
            right= right_bounds.intersect(ray, invd, hit.t, right_min, right_max);
        }
        
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
            const std::vector<unsigned int>& indices= objects[references[i].object_id].indices;
            const std::vector< const std::vector<vec3> * >& vertices= objects[references[i].object_id].vertices;            
            int slots= objects[references[i].object_id].slots;
            
            int rayt0= std::floor(ray.time * Float(slots -1));
            int rayt1= std::ceil(ray.time * Float(slots -1));
            assert(rayt0 < slots);
            assert(rayt1 < slots);
            Float t= ray.time * Float(slots -1) - Float(rayt0);
            if(t < 0) t= 0;
            if(t > 1) t= 1;
            
            int index= references[i].triangle_id;
            int v0= indices[index];
            int v1= indices[index +1];
            int v2= indices[index +2];
            
            // interpole le triangle
            Point a0= Point((*vertices[rayt0])[v0]);
            Point b0= Point((*vertices[rayt0])[v1]);
            Point c0= Point((*vertices[rayt0])[v2]);
            
            Point a1= Point((*vertices[rayt1])[v0]);
            Point b1= Point((*vertices[rayt1])[v1]);
            Point c1= Point((*vertices[rayt1])[v2]);
            
            Float rt, ru, rv;
            BVHTriangle triangle(a0 * (1 - t) + a1 * t, b0 * (1 - t) + b1 * t, c0 * (1 - t) + c1 * t, index / 3);
            if(triangle.intersect(ray, hit.t, rt, ru, rv))
            {
                hit.t= rt;   // !! raccourcir le rayon !! evite de visiter trop de noeuds apres avoir trouve une intersection
                hit.u= ru;
                hit.v= rv;
                hit.triangle_id= triangle.id;
                hit.mesh_id= references[i].object_id;
            }
        }
    }
}


bool STBVH::visible( const int index, const Ray& ray, const Vector& invd ) const
{
    assert(index != -1);
    
    const STBVHNode& node= nodes[index];
    if(ray.time < node.r0 || ray.time > node.r1) return true;

    if(node.internal())
    {
        Float t= (ray.time - node.r0) / (node.r1 - node.r0);
        if(t < 0) t= 0;
        if(t > 1) t= 1;
            
        Float tmin, tmax;
        BVHBox bounds= BVHBox(node.bounds0.pmin * (1 - t) + node.bounds1.pmin * t, node.bounds0.pmax * (1 - t) + node.bounds1.pmax * t);
        if(bounds.intersect(ray, invd, ray.tmax, tmin, tmax))
        {
            if(!visible(node.left, ray, invd))
                return false;
            if(!visible(node.right, ray, invd))
                return false;
        }
    }
    
    else // if(node.leaf())
    {
        int begin= node.leaf_begin();
        int end= node.leaf_end();
        for(int i= begin; i < end; i++)
        {
            const std::vector<unsigned int>& indices= objects[references[i].object_id].indices;
            const std::vector< const std::vector<vec3> * >& vertices= objects[references[i].object_id].vertices;            
            int slots= objects[references[i].object_id].slots;
            
            int rayt0= std::floor(ray.time * Float(slots -1));
            int rayt1= std::ceil(ray.time * Float(slots -1));
            assert(rayt0 < slots);
            assert(rayt1 < slots);
            Float t= ray.time * Float(slots -1) - Float(rayt0);
            if(t < 0) t= 0;
            if(t > 1) t= 1;
            
            int index= references[i].triangle_id;
            int v0= indices[index];
            int v1= indices[index +1];
            int v2= indices[index +2];
            
            // interpole le triangle
            Point a0= Point((*vertices[rayt0])[v0]);
            Point b0= Point((*vertices[rayt0])[v1]);
            Point c0= Point((*vertices[rayt0])[v2]);
            
            Point a1= Point((*vertices[rayt1])[v0]);
            Point b1= Point((*vertices[rayt1])[v1]);
            Point c1= Point((*vertices[rayt1])[v2]);
            
            Float rt, ru, rv;
            BVHTriangle triangle(a0 * (1 - t) + a1 * t, b0 * (1 - t) + b1 * t, c0 * (1 - t) + c1 * t, index / 3);
            if(triangle.intersect(ray, ray.tmax, rt, ru, rv))
                return false;
        }
    }
    
    return true;
}

void STBVH::intersect_notime( const int index, const Ray& ray, const Vector& invd, Hit& hit, const IntersectionFilter& filter ) const
{
    assert(index != -1);
    const STBVHNode& node= nodes[index];
    if(node.internal())
    {
        Float left_min, left_max;
        bool left= nodes[node.left].bounds0.intersect(ray, invd, hit.t, left_min, left_max);
        Float right_min, right_max;
        bool right= nodes[node.right].bounds0.intersect(ray, invd, hit.t, right_min, right_max);
        
        // parcours ordonne
        if(left && right)
        {
            Float left_mid= (left_min + left_max) / 2;
            Float right_mid= (right_min + right_max) / 2;
            if(left_mid < right_mid)
            {
                intersect_notime(node.left, ray, invd, hit, filter);
                if(hit.t >= right_min)  // ne teste l'autre fils que si necessaire...
                    intersect_notime(node.right, ray, invd, hit, filter);
            }
            else
            {
                intersect_notime(node.right, ray, invd, hit, filter);
                if(hit.t >= left_min)
                    intersect_notime(node.left, ray, invd, hit, filter);
            }
        }
        else if(left)
            intersect_notime(node.left, ray, invd, hit, filter);
        else if(right)
            intersect_notime(node.right, ray, invd, hit, filter);
    }
    
    else // if(node.leaf())
    {
        int begin= node.leaf_begin();
        int end= node.leaf_end();
        for(int i= begin; i < end; i++)
        {
            const std::vector<unsigned int>& indices= objects[references[i].object_id].indices;
            const std::vector< const std::vector<vec3> * >& vertices= objects[references[i].object_id].vertices;            
            int index= references[i].triangle_id;
            int v0= indices[index];
            int v1= indices[index +1];
            int v2= indices[index +2];
            
            // interpole le triangle
            Point a= Point((*vertices[0])[v0]);
            Point b= Point((*vertices[0])[v1]);
            Point c= Point((*vertices[0])[v2]);
            
            Float rt, ru, rv;
            BVHTriangle triangle(a, b, c, index / 3);
            if(triangle.intersect(ray, hit.t, rt, ru, rv) && filter(ray, Hit(rt, ru, rv, triangle.id, references[i].object_id)))
            {
                hit.t= rt;   // !! raccourcir le rayon !! evite de visiter trop de noeuds apres avoir trouve une intersection
                hit.u= ru;
                hit.v= rv;
                hit.triangle_id= triangle.id;
                hit.mesh_id= references[i].object_id;
            }
        }
    }
}


bool STBVH::visible_notime( const int index, const Ray& ray, const Vector& invd, const IntersectionFilter& filter  ) const
{
    assert(index != -1);
    const STBVHNode& node= nodes[index];
    if(node.internal())
    {
        Float tmin, tmax;
        if(node.bounds0.intersect(ray, invd, ray.tmax, tmin, tmax))
        {
            if(!visible_notime(node.left, ray, invd, filter))
                return false;
            if(!visible_notime(node.right, ray, invd, filter))
                return false;
        }
    }
    
    else // if(node.leaf())
    {
        int begin= node.leaf_begin();
        int end= node.leaf_end();
        for(int i= begin; i < end; i++)
        {
            const std::vector<unsigned int>& indices= objects[references[i].object_id].indices;
            const std::vector< const std::vector<vec3> * >& vertices= objects[references[i].object_id].vertices;            
            int index= references[i].triangle_id;
            int v0= indices[index];
            int v1= indices[index +1];
            int v2= indices[index +2];
            
            // interpole le triangle
            Point a= Point((*vertices[0])[v0]);
            Point b= Point((*vertices[0])[v1]);
            Point c= Point((*vertices[0])[v2]);
            
            Float rt, ru, rv;
            BVHTriangle triangle(a, b, c, index / 3);
            if(triangle.intersect(ray, ray.tmax, rt, ru, rv) && filter(ray, Hit(rt, ru, rv, triangle.id, references[i].object_id)))
                return false;
        }
    }
    
    return true;
}
