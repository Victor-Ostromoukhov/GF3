
#include <cmath>

#include "color.h"
#include "image.h"
#include "image_io.h"
#include "image_hdr.h"
#include "image_exr.h"


void range( const Image& image, float &saturation )
{
    float gmin= FLT_MAX;
    float gmax= 0.f;
    for(int i= 0; i < int(image.size()); i++)
    {
        Color32 color= image(i);
        float g= color.r + color.g + color.b;
        
        if(g < gmin) gmin= g;
        if(g > gmax) gmax= g;
    }
    
    int bins[100] = {};
    for(int i= 0; i < int(image.size()); i++)
   {
        Color32 color= image(i);
        float g= color.r + color.g + color.b;
        
        int b= (g - gmin) * 100 / (gmax - gmin);
        if(b >= 99) b= 99;
        if(b < 0) b= 0;
        bins[b]++;
    }
    
    saturation= 0;
    float qbins= 0;
    for(int i= 0; i < 100; i++)
    {
        if(qbins > .75f)
        {
            saturation= gmin + float(i) / 100 * (gmax - gmin);
            break;
        }
        
        qbins= qbins + float(bins[i]) / float(image.size());
    }
    
    printf("range [%f..%f], 75%% %f\n", gmin, gmax, saturation);
}
    

/* data/shaders/tonemap.glsl
    
    const vec3 rgby= vec3(0.3, 0.59, 0.11);
    float k1= 1.0 / pow(saturation, 1.0 / compression); // normalisation : saturation == blanc
    vec4 color= texture(image, zoom_texcoord);
    float y= dot(color.rgb, rgby);  // normalisation de la couleur : (color / y) == teinte
    color= (color / y) * k1 * pow(y, 1.0 / compression);
    
    ou 
    color.rgb= k1 * pow(color.rgb, vec3(1.0 / compression)); // sature moins les primaires... plus doux
*/
Image tone( const Image& image, const float saturation, const float gamma )
{
    Image out(image.width(), image.height());
    
    float k1= 1.f / std::pow(saturation, 1.f / gamma);
    for(int i= 0; i < int(image.size()); i++)
    {
        Color32 color= image(i);
        // float y= color.r * .3f + color.g * .59f + color.g * .11f;
        // float y= (color.r + color.g + color.g) / 3;
        // color= (color / y) * k1 * std::pow(y, 1 / gamma);
        color= Color32(std::pow(color.r, 1 / gamma), std::pow(color.g, 1 / gamma), std::pow(color.b, 1 / gamma)) * k1;
        
        out(i)= Color32(color, 1);
    }
    
    return out;
}

Image normalize_tone( const Image& image, const float saturation, const float gamma )
{
    Image out(image.width(), image.height());
    
    //~ float k1= 1.f / std::pow(saturation, 1.f / gamma);
    //~ float k1= std::pow(saturation, 1.f / gamma);
    float k1= 1 / saturation;
    for(int i= 0; i < int(image.size()); i++)
        out(i)= Color32(image(i) * k1, 1);
    
    return out;
}

Image read( const char *filename )
{
    Image image;
    if(is_exr_image(filename))
        image= read_image_exr(filename);
    else if(is_pfm_image(filename))
        image= read_image_pfm(filename);
    else if(is_hdr_image(filename))
        image= read_image_hdr(filename);
    else
        image= read_image(filename);
    
    return image;
}

int main( int argc, char *argv[] )
{
    if(argc < 2)
    {
        printf("usage: %s [--range ref.[pfm|exr|hdr]] [--exr] image.[pfm|hdr|exr]...\n", argv[0]);
        return 0;
    }
    
    float saturation= 0;
    float gamma= 2.2f;
    if(argc == 2)
    {
        Image image= read(argv[1]);
        range(image, saturation);
        write_image(tone(image, saturation, gamma), change_extension(argv[1], ".png").c_str());
    }
    else
    {
        int i= 1;
        // cherche l'option --range
        if(std::string(argv[1]) == "--range" && argc > 2)
        {
            Image ref= read(argv[2]);
            range(ref, saturation);
            i= 3;
        }
        else
        {
            // sinon, utilise la derniere image comme reference
            Image ref= read(argv[argc -1]);
            range(ref, saturation);
        }
        
        // cherche l'option --exr
        bool export_exr= false;
        if(std::string(argv[i]) == "--exr")
        {
            export_exr= true;
            i= i+1;
        }
        
        // conversion
    #pragma omp parallel for
        for(int k= i; k < argc; k++)
        {
            Image image= read(argv[k]);
            
            if(export_exr)
                write_image_exr( normalize_tone(image, saturation, gamma), change_extension(argv[k], ".exr").c_str() );
            else
                write_image( tone(image, saturation, gamma), change_extension(argv[k], ".png").c_str() );
        }
    }
    
    return 0;
}

