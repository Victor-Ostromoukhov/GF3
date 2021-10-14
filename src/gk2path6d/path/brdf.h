
#ifndef _BRDF_H
#define _BRDF_H

#include "real.h"
#include "vec.h"
#include "color.h"
#include "sampler.h"

#include "image.h"
#include "image_io.h"
#include "image_hdr.h"

#include "directions.h"

#include "options.h"


// interface pour la generation de directions sur l'hemisphere
struct Direction
{
    Direction() {}
    virtual ~Direction() {}
    
    virtual Vector sample( const Float u1, const Float u2 ) const = 0;
    virtual Float eval( const Vector& v ) const = 0;
};


// distribution 1 / 2pi
// cf GI compendium, eq 34
// http://www.cs.kuleuven.ac.be/~phil/GI/TotalCompendium.pdf
struct UniformDirection : public Direction
{
    UniformDirection( ) : Direction() {}
    
    Vector sample( const Float u1, const Float u2 ) const
    {
        Vector v;
        uniform_direction(u1, u2, v.x, v.y, v.z);
        return v;
    }
    
    // direction locale / repere de l'hemisphere    
    Float eval( const Vector& v ) const
    {
        return (v.z <= 0) ? 0 : Float(.5) / Float(M_PI);
    }
};

// distribution cos(theta) / pi
// cf GI compendium, eq 35
// http://www.cs.kuleuven.ac.be/~phil/GI/TotalCompendium.pdf
struct CosineDirection : public Direction
{
    CosineDirection( ) : Direction() {}
    
    Vector sample( const Float u1, const Float u2 ) const
    {
        Vector v;
        cosine_direction(u1, u2, v.x, v.y, v.z);
        return v;
    }
    
    // direction locale / repere de l'hemisphere    
    Float eval( const Vector& v ) const
    {
        return (v.z <= 0) ? 0 : v.z / Float(M_PI);
    }
};

//
struct BlinnDirection : public Direction
{
    BlinnDirection( const Float _alpha ) : Direction(), alpha(_alpha) {}
    
    Vector sample( const Float u1, const Float u2 ) const
    {
        Vector v;
        blinn_direction(u1, u2, alpha, v.x, v.y, v.z);
        return v;
    }
    
    Float eval( const Vector& v ) const
    {
        return (alpha + 1) / (2 * Float(M_PI)) * std::pow(v.z, alpha);        
    }
    
protected:
    Float alpha;
};

// cosine direction pdf(w)= cos theta / pi, from concentric disk mapping
// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
struct ConcentricDiskCosineDirection : public Direction
{
    ConcentricDiskCosineDirection( ) : Direction() {}
    
    Vector sample( const Float u1, const Float u2 ) const
    {
        Vector v;
        disk_cosine_direction(u1, u2, v.x, v.y, v.z);
        return v;
    }
    
    Float eval( const Vector& v ) const
    {
        return (v.z <= 0) ? 0 : v.z / Float(M_PI);
    }
};

// uniform direction pdf(w) = 1 / 2pi, from concentric disk mapping
// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
struct ConcentricDiskUniformDirection : public Direction
{
    ConcentricDiskUniformDirection( ) : Direction() {}
    
    Vector sample( const Float u1, const Float u2 ) const
    {
        Vector v;
        disk_uniform_direction(u1, u2, v.x, v.y, v.z);
        return v;
    }
    
    Float eval( const Vector& v ) const
    {
        return (v.z <= 0) ? 0 : Float(.5) / Float(M_PI);
    }
};


// blinn-phong direction pdf(w) = , from concentric disk mapping
// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
struct ConcentricDiskBlinnDirection
{
    ConcentricDiskBlinnDirection( const Float _alpha ) : alpha(_alpha) {}
    
    Vector sample( const Float u1, const Float u2 ) const
    {
        Vector v;
        disk_blinn_direction(u1, u2, alpha, v.x, v.y, v.z);
        return v;
    }

    Float eval( const Vector& v ) const
    {
        return (alpha + 1) / (2 * Float(M_PI)) * std::pow(v.z, alpha);
    }
    
protected:
    Float alpha;
};


// lambert map
struct MapDirection
{
    Vector sample( const Float u1, const Float u2 ) const
    {
        Float phi= 2 * Float(M_PI) * u1;
        Float cos_theta= 1 - u2;
        Float sin_theta= std::sqrt(std::max(Float(0), 1 - cos_theta*cos_theta));
        
        return Vector(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
    }
    
    Float eval( const Vector& v ) const
    {
        return (v.z <= 0) ? 0 : Float(.5) / Float(M_PI);
    }
};


template< typename T >
void plot_direction( const T& direction, const char *filename, const int n= 16 )
{
    FILE *out= fopen(filename, "wt");
    if(!out) return;
    
    for(int i= 0; i < n; i++)
    for(int j= 0; j < n; j++)
    {
        Float u1= Float(i + .5f) / Float(n);
        Float u2= Float(j + .5f) / Float(n);
        
        Vector w= direction.sample(u1, u2);
        fprintf(out, "%lf %lf %lf %lf %lf\n", double(u1), double(u2), double(w.x), double(w.y), double(w.z));
    }
    
    fclose(out);
}

// construit une sequence de fibonacci sur une hemisphere
// cf "Spherical Fibonacci Mapping" 
// http://lgdv.cs.fau.de/uploads/publications/spherical_fibonacci_mapping.pdf
struct Fibonacci
{
    Fibonacci( const int _n, const unsigned int _seed= 0 ) : jitter(_seed * Float(2.3283064365386963e-10)), n(_n) {}
    
    Vector direction( const int i ) const
    {
        const Float ratio= (std::sqrt(Float(5)) + 1) / 2;
        
        Float phi= 2 * Float(M_PI) * fract((i + jitter) / ratio);
        Float cos_theta= 1 - Float(2*i +1) / Float(n * 2);
        Float sin_theta= std::sqrt(1 - cos_theta*cos_theta);
        
        return Vector(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
    }
    
protected:
    Float fract( const Float v ) const
    {
        return v - std::floor(v);
    }
    
    Float jitter;
    int n;
};


// construit un repere ortho tbn, a partir d'un seul vecteur...
// cf "generating a consistently oriented tangent space" 
// http://people.compute.dtu.dk/jerf/papers/abstracts/onb.html
// cf "Building an Orthonormal Basis, Revisited", Pixar, 2017
// http://jcgt.org/published/0006/01/01/
struct World
{
    World( const Vector& _n ) : t(), b(), n(_n)
    {
        Float sign= std::copysign(Float(1), n.z);
        Float a= -1 / (sign + n.z);
        Float d= n.x * n.y * a;
        t= Vector(1 + sign * n.x * n.x * a, sign * d, -sign * n.x);
        b= Vector(d, sign + n.y * n.y * a, -n.y);        
    }
    
    // transforme le vecteur du repere local vers le repere du monde
    Vector operator( ) ( const Vector& local )  const { return local.x * t + local.y * b + local.z * n; }
    
    // transforme le vecteur du repere du monde vers le repere local
    Vector inverse( const Vector& global ) const { return Vector(dot(global, t), dot(global, b), dot(global, n)); }
    
    Vector t;
    Vector b;
    Vector n;
};


//! representation d'une distribution (isotrope) de micro facettes.
struct MF
{
    MF( const Float _alpha ) : alpha(_alpha) {}
    virtual ~MF() {}
    
    //! brdf
    virtual Float f( const Vector& wi, const Vector& wh, const Vector& wo ) const = 0;
    
    //! pdf (wi | wo)
    virtual Float eval( const Vector& wi, const Vector& wh, const Vector& wo ) const = 0;

    //! sample wi \propto pdf | wo
    virtual Vector sample( const Float u1, const Float u2, const Vector& wo ) const = 0;

protected:
    Float alpha;
};


//! MF distribution blinn : cos (theta)^ns, par defaut
struct Blinn : public MF
{
    Blinn( const Float _alpha ) : MF(_alpha) {}
    ~Blinn() {}
    
    Float f( const Vector& wi, const Vector& wh, const Vector& wo ) const
    {
        // conservation d'energie (alpha + 8) / 8pi cf "normalization zoo" http://www.thetenthplanet.de/archives/255
        Float cos_theta_h= wh.z;
        return (alpha + 8) / (8 * Float(M_PI)) * std::pow(cos_theta_h, alpha);
    }
    
    Float eval( const Vector& wi, const Vector& wh, const Vector& wo ) const
    {
        // pdf(h) + changement de mesure h / wi
        Float cos_theta_h= wh.z;
        return (alpha + 1) / (2 * Float(M_PI)) * std::pow(cos_theta_h, alpha) / (4 * dot(wo, wh));
    }
    
    Vector sample( const Float u1, const Float u2, const Vector& wo ) const
    {
        // sample h 
        Float phi= 2*Float(M_PI) * u1;
        Float cos_theta= std::pow(u2, 1 / (alpha +1));
        Float sin_theta= std::sqrt(1 - std::min(Float(1), cos_theta*cos_theta));

        Vector wh= Vector(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
        
        // reflect wi
        Vector wi= -wo + 2 * dot(wh, wo) * wh;
        return wi;
    }
};

//! representation d'une brdf : kd / pi diffus + ks micro facettes.
struct Brdf
{
    Brdf( const Float _kd, const Color& _diffuse, const Float _ks, MF& _mf, const Color& _specular )
        : mfd(_mf), dd(), diffuse(_diffuse), specular(_specular), kd(_kd), ks(_ks) {}
    
    Color f( const Vector& wi, const Vector& wo ) const
    {
        if(wi.z <= 0 || wo.z <= 0)
            return Black();
        
        Vector wh= normalize(wi + wo);
        if(wh.z <= 0) 
            return Black();
        
        if(dot(wi, wh) <= 0)
            return Black();
        if(dot(wo, wh) <= 0)
            return Black();
        
        Color d= diffuse / Float(M_PI);
        Color s= specular * mfd.f(wi, wh, wo);
        return kd*d + ks*s;
    }
    
    Vector sample( Sampler& sampler, const Vector& wo, Float& pdf )
    {
        if(options.brdf_2d)
        {
            Float u1= sampler.sample1();
            Float u2= sampler.sample2();
            return sample(u1, u2, wo, pdf);
        }
        else
        {
            Float u1= sampler.sample1();
            Float u2= sampler.sample2();
            Float u3= sampler.sample3();
            return sample(u1, u2, u3, wo, pdf);
        }
    }
    
    // utilise u1 pour selectionner le lobe de la brdf
    // utilise u1 renormalise + u2 pour generer la direction.
    Vector sample( const Float u1, const Float u2, const Vector& wo, Float& pdf ) const
    {
        pdf= 0;
        Vector wi;
        if(wo.z <= 0)
            return wi;
        
        Float k= kd / (kd + ks);
        if(u1 < k)
        {
            // 0 < u1 < k, pdf k
            Float u= u1 / k;            // renormalise u1
            wi= dd.sample(u, u2);
        }
        else
        {
            // k < u1 < k + (1-k) = 1, pdf 1-k
            Float u= (u1 - k) / (1 - k);                    // renormalise u1
            wi= mfd.sample(u, u2, wo);
        }
        
        pdf= eval(wi, wo);
        return wi;
    }
    
    // utilise u1 pour selectionner le lobe de la brdf
    // utilise u2 + u3 pour generer la direction.
    Vector sample( const Float u1, const Float u2, const Float u3, const Vector& wo, Float& pdf ) const
    {
        pdf= 0;
        Vector wi;
        if(wo.z <= 0)
            return wi;
        
        Float k= kd / (kd + ks);
        if(u1 < k)
            wi= dd.sample(u2, u3);
        else
            wi= mfd.sample(u2, u3, wo);
        
        pdf= eval(wi, wo);
        return wi;
    }
    
    Float eval( const Vector& wi, const Vector& wo ) const
    {
        if(wi.z <= 0 || wo.z <= 0)
            return 0;
        
        Vector wh= normalize(wi + wo);
        if(wh.z <= 0)
            return 0;
        
        Float d= dd.eval(wi);
        Float s= mfd.eval(wi, wh, wo);
        Float k= kd / (kd + ks);
        return k * d + (1 - k)  * s;
    }

    
    double test( Sampler& sampler, const Vector& wo= Vector(0, 0, 1) )
    {
        double sumf= 0;
        int nullpdf= 0;
        for(int i= 0; i < int(sampler.size()); i++)
        {
            Float pdf= 0;
            Vector wi= sample(sampler, wo, pdf);
            if(pdf > 0)
                sumf= sumf + f(wi, wo).r * wi.z / pdf;
            else
                nullpdf++;
        }
        sumf= sumf / double(sampler.size());
        printf("brdf %lf (null %lf)\n", sumf, double(nullpdf) / double(sampler.size()));
        
        nullpdf= 0;
        double sumpdf= 0;
        double sumff= 0;
        CosineDirection cosine;
        for(int i= 0; i < int(sampler.size()); i++)
        {
            Vector wi= cosine.sample(sampler.sample1(), sampler.sample2());
            Float pdf= cosine.eval(wi);
            if(pdf > 0)
            {
                sumff= sumff + f(wi, wo).r * wi.z / pdf;
                sumpdf= sumpdf + eval(wi, wo) / pdf;
            }
            else
                nullpdf++;
        }
        
        sumpdf= sumpdf / (double) sampler.size();
        sumff= sumff / (double) sampler.size();
        printf("brdf %lf pdf %lf (null %lf)\n", sumff, sumpdf, double(nullpdf) / double(sampler.size()));
        
        return sumf;
    }
    
    void plot( const char *filename, const Vector& wo= Vector(0, 0, 1) )
    {
        const int image_size= 1024;
        Image image(image_size, image_size);
        
        // camera Y up
        Vector Z= normalize(Vector(-1, 1, 1));
        Vector X= normalize(Vector(1, 1, 0));
        X= normalize(X - Z*dot(X, Z));
        Vector Y= - cross(X, Z);
        
        // loop over pixels
        for(int j= 0 ; j < image.height(); ++j)
        for(int i= 0 ; i < image.width() ; ++i)
        {
            // intersection point on the sphere
            Float x = Float(-1.1) + Float(2.2) * (i+Float(.5)) / (Float) image.width();
            Float y = Float(-1.1) + Float(2.2) * (j+Float(.5)) / (Float) image.height();
            if(x*x + y*y > 1)
            {
                image(i, image.height() -1 - j) = White();
                continue;
            }
            
            Float z= std::sqrt(1 - x*x - y*y);
            Vector wi= normalize(x*X + y*Y + z*Z);
            
            // evaluate function
            Float value= eval(wi, wo);
            image(i, image.height() -1 - j)= Color32(value, value, value, 1);
        }
        
        write_image_hdr(image, filename);
    }
    
protected:
    const MF& mfd;
    const CosineDirection dd;
    //~ const ConcentricDiskCosineDirection dd;
    
    Color diffuse;
    Color specular;
    Float kd;
    Float ks;
};

#endif
