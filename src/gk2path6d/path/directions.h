
#ifndef _DIRECTIONS_H
#define _DIRECTIONS_H

#include "real.h"


// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
static inline
void concentric_disk( const Float u1, const Float u2, Float& r, Float& phi )
{
    phi= 0;
    r= 0;
    
    Float a= 2 * u1 - 1;
    Float b= 2 * u2 - 1;
    
#if 0
    // cf http://psgraphics.blogspot.fr/2011/01/improved-code-for-concentric-map.html
    // not robust ?? test blinn distribution
    if(a*a > b*b)
    {
        r= a;
        phi= Float(M_PI) / 4 * (b / a);
    } 
    else 
    {
        r= b;
        if(b != 0)
            phi= Float(M_PI) / 2 - Float(M_PI) / 4 * (a / b);
    }
    
#else
    if (a > -b)
    {                    // region 1 or 2
        if (a > b)
        {                // region 1, also |a| > |b|
            r = a;
            phi = (Float(M_PI) / 4) * (b/a);
        }
        else
        {                // region 2, also |b| > |a|
            r = b;
            phi = (Float(M_PI) / 4) * (2 - (a/b));
        }
    }
    else
    {                    // region 3 or 4
        if (a < b)
        {                // region 3, also |a| >= |b|, a != 0
            r = -a;
            phi = (Float(M_PI) / 4) * (4 + (b/a));
        }
        else
        {                // region 4, |b| >= |a|, but a==0 and b==0 could occur.
            r = -b;
            if (b != 0)
                phi = (Float(M_PI) / 4) * (6 - (a/b));
        }
    }
#endif
}

static inline
// cosine direction pdf(w)= cos theta / pi, from concentric disk mapping
// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
void disk_cosine_direction( const Float u1, const Float u2, Float &x, Float &y, Float &z )
{
    Float phi= 0; 
    Float r= 0;
    concentric_disk(u1, u2, r, phi);
    
    Float u= r * std::cos(phi);
    Float v= r * std::sin(phi);
    
    x= u;
    y= v;
    z= std::sqrt(std::max(Float(0), 1 - r*r));
}

static inline
// uniform direction pdf(w)= 1 / 2pi, from concentric disk mapping
// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
void disk_uniform_direction( const Float u1, const Float u2, Float& x, Float& y, Float& z )
{
    Float phi= 0; 
    Float r= 0;
    concentric_disk(u1, u2, r, phi);
    
    Float u= r * std::cos(phi);
    Float v= r * std::sin(phi);
    
    z= std::max(Float(0), 1 - r*r);
    x= u * std::sqrt(2 - r*r);
    y= v * std::sqrt(2 - r*r);
}

static inline
// blinn direction pdf(w)= (a+1) / 2pi * cos(theta)^a, from concentric disk mapping
// "A Low Distortion Map Between Disk and Square"
// https://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf
void disk_blinn_direction( const Float u1, const Float u2, const Float alpha, Float&x, Float& y, Float& z )
{
    Float phi= 0;
    Float r= 0;
    concentric_disk(u1, u2, r, phi);
    
    Float u= r * std::cos(phi);
    Float v= r * std::sin(phi);
    
    z= std::pow(std::max(Float(0), 1 - r*r), 1 / (alpha + 1));
    x= u * std::sqrt(1 - z*z) / r;
    y= v * std::sqrt(1 - z*z) / r;
}


static inline
// uniform direction pdf(w)= 1 / 2pi
// cf GI compendium, eq 34
// http://www.cs.kuleuven.ac.be/~phil/GI/TotalCompendium.pdf
// cf pbrt3 pp 774
void uniform_direction( const Float u1, const Float u2, Float& x, Float& y, Float& z )
{
    Float phi= 2 * Float(M_PI) * u2;
    Float cos_theta= u1;
    Float sin_theta= std::sqrt(std::max(Float(0), 1 - cos_theta*cos_theta));
    
    x= std::cos(phi) * sin_theta;
    y= std::sin(phi) * sin_theta;
    z= cos_theta;
}

static inline
// cosine direction pdf(w)= cos(theta) / pi
// cf GI compendium, eq 35
// http://www.cs.kuleuven.ac.be/~phil/GI/TotalCompendium.pdf
// cf pbrt3 pp 779
void cosine_direction( const Float u1, const Float u2, Float& x, Float& y, Float& z )
{
    Float phi= 2 * Float(M_PI) * u2;
    Float cos_theta= std::sqrt(u1);
    Float sin_theta= std::sqrt(std::max(Float(0), 1 - cos_theta*cos_theta));
    
    x= std::cos(phi) * sin_theta;
    y= std::sin(phi) * sin_theta;
    z= cos_theta;
}

static inline
// blinn direction pdf(w)= pdf(w)= (a+1) / 2pi * cos(theta)^a, 
// cf GI compendium, eq 36
// http://www.cs.kuleuven.ac.be/~phil/GI/TotalCompendium.pdf
void blinn_direction( const Float u1, const Float u2, const Float alpha, Float& x, Float& y, Float& z )
{
    Float phi= 2 * Float(M_PI) * u1;
    Float cos_theta= std::pow(u2, 1 / (alpha +1));
    Float sin_theta= std::sqrt(std::max(Float(0), 1 - cos_theta*cos_theta));
    
    x= std::cos(phi) * sin_theta;
    y= std::sin(phi) * sin_theta;
    z= cos_theta;
}

#endif
