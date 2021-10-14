
#include <cstdio>
#include <algorithm>

#include "orbiter.h"

void Orbiter::lookat( const Point& center, const Float size )
{
    m_center= center;
    m_position= vec2(0, 0);
    m_rotation= vec2(0, 0);
    m_size= size;
    m_radius= size;
}

void Orbiter::lookat( const Point& pmin, const Point& pmax )
{
    lookat(center(pmin, pmax), distance(pmin, pmax));
}

void Orbiter::rotation( const Float x, const Float y )
{
    m_rotation.x= m_rotation.x + y;
    m_rotation.y= m_rotation.y + x;
}

void Orbiter::translation( const Float x, const Float y )
{
    m_position.x= m_position.x - m_size * x;
    m_position.y= m_position.y + m_size * y;
}

void Orbiter::move( const Float z )
{
    m_size= m_size - m_size * Float(0.01) * z;
    if(m_size < Float(0.01))
        m_size= Float(0.01);
}

Transform Orbiter::view( ) const
{
    return Translation( -m_position.x, -m_position.y, -m_size ) 
        * RotationX(m_rotation.x) * RotationY(m_rotation.y) 
        * Translation( - Vector(m_center) ); 
}

Transform Orbiter::projection( const Float width, const Float height, const Float fov ) const
{
    // calcule la distance entre le centre de l'objet et la camera
    //~ Transform t= view();
    //~ Point c= t(m_center);
    //~ Float d= -c.z;
    Float d= distance(m_center, Point(m_position.x, m_position.y, m_size));     // meme resultat plus rapide a calculer
    
    // regle near et far en fonction du centre et du rayon englobant l'objet 
    return Perspective(fov, width / height, std::max(Float(0.1), d - m_radius), std::max(Float(1), d + m_radius));
}

void Orbiter::frame( const Float width, const Float height, const Float z, const Float fov, Point& dO, Vector& dx, Vector& dy ) const
{
    Transform v= view();
    Transform p= projection(width, height, fov);
    Transform viewport= Viewport(width, height);
    Transform t= viewport * p * v;              // passage monde vers image
    Transform tinv= t.inverse();                // l'inverse, passage image vers monde
    
    // origine du plan image
    dO= tinv(Point(0, 0, z));
    // axe x du plan image
    Point d1= tinv(Point(1, 0, z));
    // axe y du plan image
    Point d2= tinv(Point(0, 1, z));
    
    dx= Vector(dO, d1);
    dy= Vector(dO, d2);
}

Point Orbiter::position( )
{
    Transform t= view();     // passage monde vers camera
    Transform tinv= t.inverse();            // l'inverse, passage camera vers monde
    
    return tinv(Point(0, 0, 0));        // la camera se trouve a l'origine, dans le repere camera...
}

int Orbiter::read_orbiter( const char *filename )
{
    FILE *in= fopen(filename, "rt");
    if(in == nullptr)
    {
        printf("[error] loading orbiter '%s'...\n", filename);
        return -1;
    }
    
    printf("loading orbiter '%s'...\n", filename);
    
    bool errors= false;
    double tmp[4]= {};
    
    if(fscanf(in, "c %lf %lf %lf \n", &tmp[0], &tmp[1], &tmp[2]) != 3)
        errors= true;
    m_center= Point(tmp[0], tmp[1], tmp[2]);
    
    if(fscanf(in, "p %lf %lf \n", &tmp[0], &tmp[1]) != 2)
        errors= true;
    m_position= vec2(tmp[0], tmp[1]);
    
    if(fscanf(in, "r %lf %lf \n", &tmp[0], &tmp[1]) != 2)
        errors= true;
    m_rotation= vec2(tmp[0], tmp[1]);
    
    if(fscanf(in, "s %lf %lf \n", &tmp[0], &tmp[1]) != 2)
        errors= true;
    m_size= tmp[0];
    m_radius= tmp[1];
    
    fclose(in);
    if(errors)
    {
        printf("[error] loading orbiter '%s'...\n", filename);
        return -1;
    }
    
    return 0;
}

int Orbiter::write_orbiter( const char *filename )
{
    FILE *out= fopen(filename, "wt");
    if(out == nullptr)
        return -1;
    
    printf("writing orbiter '%s'...\n", filename);
    
    fprintf(out, "c %lf %lf %lf\n", double(m_center.x), double(m_center.y), double(m_center.z));
    fprintf(out, "p %lf %lf\n", double(m_position.x), double(m_position.y));
    fprintf(out, "r %lf %lf\n", double(m_rotation.x), double(m_rotation.y));
    fprintf(out, "s %lf %lf\n", double(m_size), double(m_radius));
    
    fclose(out);
    return 0;
}
