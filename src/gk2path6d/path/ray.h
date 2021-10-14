
#ifndef _RAY_H
#define _RAY_H

#include "vec.h"


struct Ray
{
    Point o;
    Float time;
    Vector d;
    Float tmax;
    
    Ray( ) : o(), time(0), d(), tmax(0) {}
    
    Ray( const Point origine, const Point extremite ) : o(origine), time(0), d(Vector(origine, extremite)), tmax(1) {}
    Ray( const Point origine, const Vector direction ) : o(origine), time(0), d(direction), tmax(FLOAT_MAX) {}
    
    Ray( const Float _time, const Point _origine, const Point _extremite ) : o(_origine), time(_time), d(Vector(_origine, _extremite)), tmax(1) {}
    Ray( const Float _time, const Point _origine, const Vector _direction ) : o(_origine), time(_time), d(_direction), tmax(FLOAT_MAX) {}
    
    Point operator( ) ( const Float t ) const { return o + t * d; }
};

struct Hit
{
    Float t, u, v;
    int triangle_id;
    int mesh_id;
    
    Hit( ) : t(FLOAT_MAX), u(0), v(0), triangle_id(-1), mesh_id(-1) {}
    Hit( const Float& _t, const Float _u, const Float& _v, const int _tid, const int _mid ) : t(_t), u(_u), v(_v), triangle_id(_tid), mesh_id(_mid) {}
};


struct Raydata
{
    Ray ray;
    Hit hit;
    
    Raydata( ) : ray(), hit() {}
    Raydata( const Ray& r, const Hit& h ) : ray(r), hit(h) {}
    
    Raydata( const Point& p, const Point& q ) : ray(p, q), hit() {}
    Raydata( const Point& o, const Vector& d ) : ray(o, d), hit() {}
    Raydata( const Point& p, const Point& q, const Float t ) : ray(t, p, q), hit() {}
    Raydata( const Point& o, const Vector& d, const Float t ) : ray(t, o, d), hit() {}
    
    Point o( ) const { return ray.o; }
    Vector d( ) const { return normalize(ray.d); }
    Float time( ) const { return ray.time; }
    
    bool intersect( ) const { return (hit.mesh_id != -1); }
    
    int index( ) const { return hit.mesh_id; }
    int hitid( ) const { return hit.triangle_id; }
    Point hitp( ) const { assert(intersect()); return ray(hit.t); }
    Float hitu( ) const { assert(intersect()); return hit.u; }
    Float hitv( ) const { assert(intersect()); return hit.v; }
};

#endif
