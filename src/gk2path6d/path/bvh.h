
#ifndef _BVH_H
#define _BVH_H

#include <cassert>
#include <vector>

#include "real.h"
#include "vec.h"
#include "ray.h"
#include <algorithm>

//! representation d'une boite englobante alignee sur les axes.
struct BVHBox
{
    Point pmin;
    Point pmax;
    
    BVHBox( ) : pmin(FLOAT_MAX, FLOAT_MAX, FLOAT_MAX), pmax(-FLOAT_MAX, -FLOAT_MAX, -FLOAT_MAX) {}
    BVHBox( const Point& a, const Point& b ) : pmin(min(a, b)), pmax(max(a, b)) {}
    BVHBox( const BVHBox& a, const BVHBox& b ) : pmin(min(a.pmin, b.pmin)), pmax(max(a.pmax, b.pmax)) {}
    
    void clear( ) 
    {
        pmin= Point(FLOAT_MAX, FLOAT_MAX, FLOAT_MAX);
        pmax= Point(-FLOAT_MAX, -FLOAT_MAX, -FLOAT_MAX);
    }

    // ajoute un point
    void insert( const Point& p )
    {
        pmin= min(pmin, p);
        pmax= max(pmax, p);
    }
    
    // ajoute une bbox
    void insert( const BVHBox& b )
    {
        pmin= min(pmin, b.pmin);
        pmax= max(pmax, b.pmax);
    }
    
    bool empty( ) const
    {
        bool x= (pmax.x < pmin.x);
        bool y= (pmax.y < pmin.y);
        bool z= (pmax.z < pmin.z);
        return (x || y || z);
    }
    
    // renvoie l'aire de la bbox.
    Float area( ) const
    {
        if(empty()) return 0;
        
        Vector d(pmin, pmax);
        return 2 * d.x * d.y + 2 * d.x * d.z + 2 * d.y * d.z;
    }
    
    Point center( ) const
    {
        return (pmin + pmax) / 2;
    }
    
    bool inside( const Point& p ) const
    {
        return (p.x >= pmin.x && p.x <= pmax.x
            && p.y >= pmin.y && p.y <= pmax.y
            && p.z >= pmin.z && p.z <= pmax.z);
    }
    
    bool intersect( const Ray& ray, const Vector& invd, const Float htmax, Float& rtmin, Float &rtmax ) const
    {
        Point rmin= pmin;
        Point rmax= pmax;
        if(ray.d.x < 0) std::swap(rmin.x, rmax.x);
        if(ray.d.y < 0) std::swap(rmin.y, rmax.y);
        if(ray.d.z < 0) std::swap(rmin.z, rmax.z);
        Vector dmin= (rmin - ray.o) * invd;
        Vector dmax= (rmax - ray.o) * invd;
        
        rtmin= std::max(dmin.z, std::max(dmin.y, std::max(dmin.x, Float(0))));
        rtmax= std::min(dmax.z, std::min(dmax.y, std::min(dmax.x, htmax)));
        /*  cf "Robust BVH Ray Traversal", T. Ize, 
            http://www.cs.utah.edu/~thiago/papers/robustBVH-v2.pdf */
        //~ rtmax= rtmax * 1.00000024f;
        return (rtmin <= rtmax);
    }
};

//! representation d'un triangle "geometrique".
struct BVHTriangle
{
    Point p;
    Vector e1;
    Vector e2;
    int id;
    
    BVHTriangle( const Point& _a, const Point& _b, const Point& _c, const int _id ) : p(_a), e1(Vector(_a, _b)), e2(Vector(_a, _c)), id(_id) {}
    BVHTriangle( const BVHTriangle& _t, const int _id ) : p(_t.p), e1(_t.e1), e2(_t.e2), id(_id) {}
    
    Point a( ) const { return p; }
    Point b( ) const { return p + e1; }
    Point c( ) const { return p + e2; }
    
    Vector normal( ) const { return normalize(cross(e1, e2)); }
    Point point( const Float u, const Float v ) const { return p + u * e1 + v * e2;}
    Float area( ) const { return length(cross(e1, e2)) / 2; }
    
    BVHBox bounds( ) const
    {
        BVHBox box;
        box.insert(a());
        box.insert(b());
        box.insert(c());
        return box;
    }
    
    /* calcule l'intersection ray/triangle
        cf "fast, minimum storage ray-triangle intersection" 
        http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
        
        renvoie faux s'il n'y a pas d'intersection valide, une intersection peut exister mais peut ne pas se trouver dans l'intervalle [0 htmax] du rayon. \n
        renvoie vrai + les coordonnees barycentriques (ru, rv) du point d'intersection + sa position le long du rayon (rt). \n
        convention barycentrique : p(u, v)= (1 - u - v) * a + u * b + v * c \n
    */
    bool intersect( const Ray &ray, const Float htmax, Float &rt, Float &ru, Float&rv ) const
    {
        Vector pvec= cross(ray.d, e2);
        Float det= dot(e1, pvec);
        
        Float inv_det= 1 / det;
        Vector tvec(p, ray.o);

        Float u= dot(tvec, pvec) * inv_det;
        if(u < 0 || u > 1) return false;

        Vector qvec= cross(tvec, e1);
        Float v= dot(ray.d, qvec) * inv_det;
        if(v < 0 || u + v > 1) return false;

        rt= dot(e2, qvec) * inv_det;
        ru= u;
        rv= v;
        return (rt <= htmax && rt >= 0);
    }
};


//! representation d'un noeud / feuille d'un bvh.
struct BVHNode 
{
    BVHBox bounds;
    int left;           //!< fils gauche et droit pour les noeuds internes, right < 0 pour les feuilles
    int right;          //!< feuilles referencent une sequence de triangles [left .. left - right[
    
    // constructeur nomme, construit une feuille
    BVHNode& make_leaf( const BVHBox& _box, const int _begin, const int _end )
    {
        assert(_begin < _end);
        bounds= _box;
        left= _begin;
        right= - (_end - _begin);       // - n
        assert(right < 0);
        return *this;
    }
    
    // constructeur nomme, construit un noeud interne
    BVHNode& make_node( const BVHBox& _box, const int _left, const int _right )
    {
        bounds= _box;
        left= _left;
        right= _right;
        return *this;
    }
    
    // renvoie vrai si le noeud est interne ou une feuille
    bool internal( ) const { return (right >= 0); }
    bool leaf( ) const { return (right < 0); }
    int leaf_begin( ) const { assert(leaf()); return left; }
    int leaf_n( ) const { assert(leaf()); return -right; }
    int leaf_end( ) const { assert(leaf()); return left - right; }
};


//! representation d'un triangle/objet pour construire un bvh
struct BVHPrimitive
{
    BVHBox bounds;
    Point center;
    int object_id;
    
    BVHPrimitive( const BVHBox& _box, const int _id ) : bounds(_box), center(_box.pmin + (_box.pmax - _box.pmin) / 2), object_id(_id) {}
};


//! representation d'une hierarchie d'englobants.
struct BVH
{
    std::vector<BVHNode> nodes;
    std::vector<BVHTriangle> triangles;
    std::vector<int> node_axis;
    BVHBox bounds;
    int root;
    
    BVH( );
    void build( );      //!< attention : les triangles doivent etre affectes avant !!
    
    bool intersect( const Ray& ray, Hit& hit ) const;
    bool visible( const Point& p, const Point& q ) const;

protected:
    void intersect( const int index, const Ray& ray, const Vector& invd, Hit& hit ) const;
    bool visible( const int index, const Ray& ray, const Vector& invd ) const;
};


#endif
