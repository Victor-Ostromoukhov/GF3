
#ifndef _CAMERA_H
#define _CAMERA_H

#include "real.h"
#include "directions.h"


// profondeur de champ de la camera
struct Lens
{
    Lens( ) : radius(0), focal(0) {}
    Lens( const Float _radius, const Float _focal ) : radius(_radius), focal(_focal) {}
    
    // o, d : origine et direction, repere camera, o= Origin()
    // lens_o, lens_d : origine direction, repere camera
    Float sample( Sampler& sampler, const Point& o, const Vector& d, Point& lens_o, Vector& lens_d ) const
    {
        return sample(sampler.sample1(), sampler.sample2(), o, d, lens_o, lens_d);
    }
    
    Float sample( const Float u1, const Float u2, const Point& vo, const Vector& vd, Point& lens_o, Vector& lens_d ) const
    {
        lens_o= vo;
        lens_d= vd;
        if(radius <= 0) 
            return 1;
        
    #ifndef LENS_SQUARE
        // sample disk
        Float r, phi;
        concentric_disk(u1, u2, r, phi);
        lens_o= radius * r * Point(std::cos(phi), std::sin(phi), 0);
    #else
        #warning lens: using square sampling
        // sample square
        lens_o= radius * Point(2*u1 -1, 2*u2 -1, 0);
    #endif
        
        Vector d= normalize(vd);
        Float ft= focal / -d.z;
        lens_d= Point(ft * d) - lens_o;
        
        return eval(lens_o, lens_d);
    }
    
    Float eval( const Point& vo, const Vector& vd ) const
    {
        if(radius > 0)
            return 1 / -normalize(vd).z;
        else
            return 1;
    }
    
    Float area( ) const 
    { 
    #ifndef LENS_SQUARE
        // sample disk
        return radius*radius * Float(M_PI); 
    #else
        // sample square
        return 4 * radius*radius;
    #endif
    }
    
    //! relit depuis un fichier texte. 
    int read_lens( const char *filename )
    {
        FILE *in= fopen(filename, "rt");
        if(in == nullptr)
        {
            printf("[error] loading lens '%s'...\n", filename);
            return -1;
        }
        
        printf("loading lens '%s'...\n", filename);
        
        double r, f;
        bool errors= false;
        if(fscanf(in, "l %lf %lf \n", &r, &f) != 2)
            errors= true;
        
        fclose(in);
        if(errors)
        {
            printf("[error] loading lens '%s'...\n", filename);
            return -1;
        }
        
        radius= r;
        focal= f;
        return 0;
    }
    
    //! enregistre dans un fichier texte.
    int write_lens( const char *filename )
    {
        FILE *out= fopen(filename, "wt");
        if(out == nullptr)
            return -1;
        
        printf("writing lens '%s'...\n", filename);
        
        fprintf(out, "l %lf %lf\n", double(radius), double(focal));
        
        fclose(out);
        return 0;
    }
    
protected:
    Float radius;
    Float focal;
};


struct Camera
{
    Camera( ) : 
        view(), projection(), viewport(),
        i2c(), c2w(), lens(), fov(0)
    {}
    
    Camera( const int _width, const int _height, const Float _fov, const Orbiter& _orbiter, const Lens& _lens ) 
        : 
        view(_orbiter.view()), 
        projection(_orbiter.projection(_width, _height, _fov)), 
        viewport(Viewport(_width, _height)),
        i2c(Inverse(viewport * projection)),
        c2w(Inverse(view)),
        lens(_lens),
        fov(_fov)
        //~ , pixel_cone_angle(0)
    {
        Float fovx= 2 * radians(fov);
        Float fovy= Float(_width) / Float(_height) * fovx;
        //~ pixel_cone_angle= std::atan(2 * std::tan(fovy / 2) / Float(_height));
    }
    
    // image sample (u1, u2) in [0 width]x[0 height], (u3, u4) in [0 1]x[0 1]
    Float sample( Sampler& sampler, Point& o, Vector& d ) 
    { 
        return sample(sampler.sample1(), sampler.sample2(), sampler.sample3(), sampler.sample4(), o, d); 
    }
    
    // (u1, u2) in [0 width] x [0 height]
    Float sample( const Float u1, const Float u2, const Float u3, const Float u4, Point& o, Vector& d )
    {
        // repere image, plan far
        Point ie= Point(u1, u2, 1);  
        
        // passage repere camera
        Point ce= i2c(ie);
        Point co= Origin();
        Vector cd(co, ce);
        
        // perturbation par l'optique
        Point lo;
        Vector ld;
        Float pdf= lens.sample(u3, u4, co, cd, lo, ld);
        
        // passage repere monde
        o= c2w(lo);
        d= c2w(ld);
        return pdf;
    }
    
    Float eval( const Point& vo, const Vector& vd ) { return lens.eval(vo, vd); }

    //~ Float pixel_spread( ) const { return pixel_cone_angle; }
    
protected:
    Transform view;
    Transform projection;
    Transform viewport;
    Transform i2c;      // image to camera
    Transform c2w;      // camera to world
    Lens lens;
    Float fov;
    //~ Float pixel_cone_angle;
};

#endif
