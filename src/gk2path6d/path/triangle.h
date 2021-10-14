
#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include "vec.h"
#include "mesh.h"


struct Triangle
{
    Point a, b, c;
    Vector na, nb, nc;
    vec2 ta, tb, tc;
    int ida, idb, idc;
    
    Triangle( ) = default;
    
    Triangle( const TriangleData &data ) 
        : 
        a(Point(data.a)), b(Point(data.b)), c(Point(data.c)),
        na(Vector(data.na)), nb(Vector(data.nb)), nc(Vector(data.nc)),
        ta(data.ta), tb(data.tb), tc(data.tc),
        ida(data.ida), idb(data.idb), idc(data.idc)
    {}
    
    //! renvoie l'aire du triangle
    Float area( ) const
    {
        return length(cross(Vector(a, b), Vector(a, c))) / 2;
    }
    
    //! renvoie la normale du triangle
    Vector normal( ) const
    {
        return normalize(cross(Vector(a, b), Vector(a, c)));
    }
    
    //! renvoie un point a l'interieur du triangle connaissant ses coordonnees barycentriques.
    //! convention p(u, v)= (1 - u - v) * a + u * b + v * c
    Point point( const Float u, const Float v ) const
    {
        Float w= 1 - u - v;
        return a * w + b * u + c * v;
    }
    
    //! renvoie les texcoord, connaissant les coordonnees barycentriques. 
    //! convention p(u, v)= (1 - u - v) * a + u * b + v * c
    vec2 texcoord( const Float u, const Float v ) const
    {
        Float w= 1 - u - v;
        return vec2(ta.x * w + tb.x * u + tc.x * v, ta.y * w + tb.y * u + tc.y * v);
    }
    
    //! renvoie une normale a l'interieur du triangle connaissant ses coordonnees barycentriques.
    //! convention p(u, v)= (1 - u - v) * a + u * b + v * c
    Vector normal( const Float u, const Float v ) const
    {
        Float w= 1 - u - v;
        return normalize(na * w + nb * u + nc * v);
    }
};

#endif
