
#ifndef _POINTS_H
#define _POINTS_H

#include <cmath>

#include "real.h"
#include "vec.h"
#include "color.h"
#include "triangle.h"
#include "sampler.h"

#include "options.h"


struct Source
{
    Color emission;
    Triangle triangle;
    
    Source( ) : emission(), triangle() {}
    Source( const TriangleData& data, const Color& color ) : emission(color), triangle(data) {}
    
    Float area( ) const { return triangle.area(); }
    Vector normal( ) const { return triangle.normal(); }
    
    Point sample( const Float u1, const Float u2 ) const
    {
    #ifndef TRIANGLE_MAPPING
#ifdef _MSC_VER
        #pragma message ("Warning : points.h  sqrt mapping" )
#else
        #warning source: sqrt mapping
#endif
        // GI compendium eq 18
        Float r1= std::sqrt(u1);
        Float r2= u2;
        
        Float b0= (1 - r1);
        Float b1= (1 - r2) * r1;
        Float b2= r2 * r1;
        
    #else
        // A Low-Distortion Map Between Triangle and Square
        // cf https://pharr.org/matt/blog/2019/03/13/triangle-sampling-1.5.html
        // cf https://drive.google.com/file/d/1J-183vt4BrN9wmqItECIjjLIKwm29qSg/view
        Float b0= u1 / 2;
        Float b1= u2 / 2;
        Float offset= b1 - b0;
        if(offset > 0)
            b1= b1 + offset;
        else
            b0= b0 - offset;
        Float b2= 1 - b0 - b1;
    #endif
        
        // select smallest area vertex
        Point a= triangle.a;
        Point b= triangle.b;
        Point c= triangle.c;
        
        Float fa= distance2(a, b) * distance2(a, c);
        Float fb= distance2(b, a) * distance2(b, c);
        if(fb < fa)
        {
            a= triangle.b;
            b= triangle.c;
            c= triangle.a;
        }
        
        Float fc= distance2(c, a) * distance2(c, b);
        if(fc < fb && fc < fa)
        {
            a= triangle.c;
            b= triangle.a;
            c= triangle.b;
        }
        
        Point s= a * b0 + b * b1 + c * b2;
        return s;
    }
};

struct QuadSource
{
    Color emission;
    Point a, b, c, d;
    
    QuadSource( ) : emission(), a(), b(), c(), d() {}
    QuadSource( const Point& _a, const Point& _b, const Point& _c, const Point& _d, const Color& color ) : emission(color), a(_a), b(_b), c(_c), d(_d) {}
    
    Float area( ) const { return length(cross(Vector(a, b), Vector(a, d))); }
    Vector normal( ) const { return normalize(cross(Vector(a, b), Vector(a, d))); }
    
    Point sample( const Float u1, const Float u2 ) const
    {
        Point s= a * ((1 - u1) * (1 - u2))
               + b * (u1       * (1 - u2))
               + d * ((1 - u1) * u2)
               + c * (u1       * u2);
        return s;
    }    
};


struct LightSample
{
    Point s;
    Vector sn;
    Color emission;
    int source_id;
    
    LightSample( ) : s(), sn(), emission(), source_id(-1) {}
    LightSample( const Point& _s, const Vector& _sn, const Color& _emission, const int _id= -1 ) : s(_s), sn(_sn), emission(_emission), source_id(_id) {}
};

static
Vector normal( const Point& a, const Point& b, const Point& c )
{
    return normalize(cross(Vector(a, b), Vector(a, c)));
}


// choisit une source 1/A.
struct UniformSource
{
    std::vector<Source> sources;
    std::vector<QuadSource> quads;
    
    UniformSource( ) = default;
    
    UniformSource( const std::vector<Source>& _sources ) : sources(_sources), quads(), cdf(), total_area(0), total_emission(0)
    {
        if(options.source_quads || options.source_quad1)
        {
            // reconstruit les quads, si possible
            for(int i= 0; i +1< int(sources.size()); i++)
            {
                for(int e= 0; e < 3; e++)
                {
                    // origine et extremite de l'arete e dans la source
                    int e1= (&sources[i].triangle.ida)[e];
                    int e2= (&sources[i].triangle.ida)[(e+1) %3];
                    if(e2 < e1) std::swap(e1, e2);  // ordonne les aretes
                    
                    // recherche l'arete dans la source suivante
                    for(int k= 0; k < 3; k++)
                    {
                        int k1= (&sources[i+1].triangle.ida)[k];
                        int k2= (&sources[i+1].triangle.ida)[(k+1) %3];
                        if(k2 < k1) std::swap(k1, k2);      // ordonne les aretes
                        
                        if(e1 == k1 && e2 == k2)
                        {
                            // reconstruit les indices des aretes
                            int e1= (&sources[i].triangle.ida)[e];
                            int e2= (&sources[i].triangle.ida)[(e+1) %3];
                            int k1= (&sources[i+1].triangle.ida)[k];
                            int k2= (&sources[i+1].triangle.ida)[(k+1) %3];
                            int e3= (&sources[i].triangle.ida)[(e+2) %3];
                            int k3= (&sources[i+1].triangle.ida)[(k+2) %3];
                            //~ printf("quad %d %d: e %d %d %d, k %d %d %d\n", i, i+1, e, e1, e2, k, k1, k2);
                            
                            if(e2 == sources[i].triangle.ida)
                            {
                                //~ printf("  e2= a: %d %d %d %d, areas %f %f\n", e2, e3, e1, k3, 
                                    //~ dot(normal((&sources[i].triangle.a)[(e+1)%3], (&sources[i].triangle.a)[(e+2)%3], (&sources[i].triangle.a)[e]), sources[i].normal()),
                                    //~ dot(normal((&sources[i].triangle.a)[(e+1)%3], (&sources[i].triangle.a)[e], (&sources[i+1].triangle.a)[(k+2)%3]), sources[i].normal()));
                                
                                quads.push_back( QuadSource(
                                        (&sources[i].triangle.a)[(e+1)%3], 
                                        (&sources[i].triangle.a)[(e+2)%3], 
                                        (&sources[i].triangle.a)[e], 
                                        (&sources[i+1].triangle.a)[(k+2)%3], sources[i].emission) );
                            }                        
                            else if(e2 == sources[i].triangle.idb)
                            {
                                //~ printf("  e2= b: %d %d %d %d, areas %f %f\n", e1, k3, e2, e3,
                                    //~ dot(normal((&sources[i].triangle.a)[e], (&sources[i+1].triangle.a)[(k+2)%3], (&sources[i].triangle.a)[(e+1)%3]), sources[i].normal()),
                                    //~ dot(normal((&sources[i].triangle.a)[e], (&sources[i].triangle.a)[(e+1)%3], (&sources[i].triangle.a)[(e+2)%3]), sources[i].normal()));
                                
                                quads.push_back( QuadSource(
                                    (&sources[i].triangle.a)[e], 
                                    (&sources[i+1].triangle.a)[(k+2)%3], 
                                    (&sources[i].triangle.a)[(e+1)%3], 
                                    (&sources[i].triangle.a)[(e+2)%3], sources[i].emission) );
                            }
                            
                            //~ else
                                //~ printf("  oops\n");
                        }
                    }
                }
            }
            
            assert(quads.size() >= 1);
            if(options.source_quad1)
                // force l'utilisation d'une seule source
                quads.resize(1);
            
            total_area= 0;
            total_emission= 0;
            for(int i= 0; i < int(quads.size()); i++)
            {
                Float v= quads[i].area();
                assert(v > 0);
                
                total_area= total_area + v;
                total_emission= total_emission + v * quads[i].emission.power();
                cdf.push_back(total_area);
                
                //~ printf("  quad source %d: area %f, emission %f\n", i, sources[i].area(), sources[i].emission.grey() * sources[i].area());
            }
            
            printf("%d quad sources, %f area, %f emission\n", size(), total_area, total_emission);
        }
        else
        {
            total_area= 0;
            total_emission= 0;
            for(int i= 0; i < int(sources.size()); i++)
            {
                Float v= sources[i].area();
                assert(v > 0);
                
                total_area= total_area + v;
                total_emission= total_emission + v * sources[i].emission.power();
                cdf.push_back(total_area);
                
                //~ printf("  source %d: area %f, emission %f\n", i, sources[i].area(), sources[i].emission.grey() * sources[i].area());
            }
            
            printf("%d sources, %f area, %f emission\n", size(), total_area, total_emission);
        }
    }
    
    int size( ) const
    {
        if(options.source_quads || options.source_quad1)
            return int(quads.size());
        else
            return int(sources.size());
    }
    
    LightSample sample( Sampler& sampler, Float& pdf ) const
    {
        if(options.source_2d)
        {
            Float u1= sampler.sample1();
            Float u2= sampler.sample2();
            return sample(u1, u2, pdf);
        }
        else
        {
            Float u1= sampler.sample1();
            Float u2= sampler.sample2();
            Float u3= sampler.sample3();
            return sample(u1, u2, u3, pdf);
        }
    }
    
    LightSample sample( const Float u1, const Float u2, Float& pdf ) const
    {
        int id= -1;
        // choisit une source u1
        Float r= u1 * total_area;

        // recherche dichotomique
        int p= 0;
        int q= int(cdf.size()) -1;
        while(p < q)
        {
            int m= (p+q) / 2;
            if(cdf[m] < r)
               p= m+1;
            else
               q= m;
        }
        assert(p >= 0 && p < (int) cdf.size());
        id= p;
        
        // renormalise u1 : cdf[p-1] < r < cdf[p]
        Float pp= 0;
        if(p > 0) pp= cdf[p-1];
        Float u= (r - pp) / (cdf[p] - pp);
        assert(u >= 0 && u <= 1);
        
        // choisit un point sur la source (u, u2)
        pdf= 1 / total_area;
        if(options.source_quads || options.source_quad1)
            return LightSample(quads[id].sample(u, u2), quads[id].normal(), quads[id].emission, id);
        else
            return LightSample(sources[id].sample(u, u2), sources[id].normal(), sources[id].emission, id);
    }

    LightSample sample( const Float u1, const Float u2, const Float u3, Float& pdf ) const
    {
        int id= -1;
        // choisit une source en fonction de son aire
        Float r= u1 * total_area;

        // recherche dichotomique
        int p= 0;
        int q= int(cdf.size()) -1;
        while(p < q)
        {
            int m= (p+q) / 2;
            if(cdf[m] < r)
               p= m+1;
            else
               q= m;
        }
        assert(p >= 0 && p < int(cdf.size()));
        id= p;
        
        // choisit un point sur la source (u2, u3)
        pdf= 1 / total_area;
        if(options.source_quads || options.source_quad1)
            return LightSample(quads[id].sample(u2, u3), quads[id].normal(), quads[id].emission, id);
        else
            return LightSample(sources[id].sample(u2, u3), sources[id].normal(), sources[id].emission, id);
    }
    
    Float eval( const Point& p, const Material& material ) const
    {
        // verifier que le point se trouve sur une source
        if(material.emission.power() == 0) 
            return 0;
        
        return 1 / total_area;
    }
    
protected:
    std::vector<Float> cdf;
    Float total_area;
    Float total_emission;
};

#endif
