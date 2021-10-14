
#include <omp.h>

#include <cstdio>
#include <cstring>
#include <utility>
#include <algorithm>

#include "image.h"
#include "image_io.h"
#include "image_hdr.h"
#include "image_exr.h"
#include "color.h"

#include "llsq.hpp"

#include "colormaps/colormapTemperatureMap.h"
const unsigned char *colormap8ui= colormapTemperatureMap;

Color32 colormap( const double x, const double xmax )
{
    int id= 512 * x / xmax;
    if(id < 0) id= 0;
    if(id >= 512) id= 511;
    
    return Color32(colormap8ui[id*3 +0] / 255.0, colormap8ui[id*3 +1] / 255.0, colormap8ui[id*3 +2] / 255.0);
}

//~ double MAXMSE = 25.0; // simple
//~ double MAXMSE = 0.002;     // cornell
double MAXMSE = 0.001;     // cornell indirect
double MAXOFFSET = 200.0;

#if 0
std::pair<double,double> MSE( const Image &A, const Image &B, double &vmin, double &vmax )
{
    double val= 0.0;
    double var= 0.0;
    vmax= -1.0;
    vmin= 100000.0;
    
    for(size_t i= 0; i < A.size(); ++i)
    {
        const Color32& a= A(i);
        const Color32& b= B(i);
        
        Color32 d= A(i) - B(i);
        double v= (d.r*d.r + d.g*d.g + d.b*d.b) / 3;
        val += v;
        vmax= std::max(vmax, v);
        vmin= std::min(vmin, v);
        var += v*v;
    }
    
    val /= double(A.size() * A.size());
    var= (var - val*val);
    return std::make_pair(val, var) ;
}
#endif

std::pair<double, double> MSElocal(const Image &A, const Image &B,
    int mX, int mY,
    int MX, int MY ,
    double &vmin, double &vmax)
{
    double val= 0.0;
    double var= 0.0;
    vmax= -1.0;
    vmin= 100000.0;
    
    int cpt= 0;
    for(int i=mX; i <= MX; ++i)
    for(int j=mY; j <= MY; ++j)
    {
        Color32 d= A(i, j) - B(i, j);
        double v= (d.r*d.r + d.g*d.g + d.b*d.b) / 3;
        val += v;
        vmax= std::max(vmax, v );
        vmin= std::min(vmin, v );
        var += v*v;
        cpt++;
    }
    
    val /=double(cpt);
    var /= double(cpt);
    var= (var - val*val);
    return std::make_pair(val, var);
}

std::pair<double, double> Mlocal( /* ref */ const Image &A, /* image */ const Image &B,
    int mX, int mY,
    int MX, int MY ,
    double &vmin, double &vmax)
{
    double val= 0.0;
    double var= 0.0;
    vmax= -1.0;
    vmin= 100000.0;
    
    int cpt= 0;
    for(int i=mX; i <= MX; ++i)
    for(int j=mY; j <= MY; ++j)
    {
        Color32 d= A(i, j) - B(i, j);
        //~ double v= (d.r*d.r + d.g*d.g + d.b*d.b) / 3;
        double v= (std::abs(d.r) + std::abs(d.g) + std::abs(d.b)) / 3;
        val += v;
        vmax= std::max(vmax, v );
        vmin= std::min(vmin, v );
        var += v*v;
        cpt++;
    }
    
    val /=double(cpt);
    var /= double(cpt);
    var= (var - val*val);
    return std::make_pair(val, var);
}

std::pair<double, double> Mlocal_error( /* image */ const Image &B,
    int mX, int mY,
    int MX, int MY ,
    double &vmin, double &vmax)
{
    double val= 0.0;
    double var= 0.0;
    vmax= -1.0;
    vmin= 100000.0;
    
    int cpt= 0;
    for(int i=mX; i <= MX; ++i)
    for(int j=mY; j <= MY; ++j)
    {
        //~ Color32 d= A(i, j) - B(i, j);
        Color32 d= B(i, j);
        //~ double v= (d.r*d.r + d.g*d.g + d.b*d.b) / 3;
        double v= (std::abs(d.r) + std::abs(d.g) + std::abs(d.b)) / 3;
        val += v;
        vmax= std::max(vmax, v );
        vmin= std::min(vmin, v );
        var += v*v;
        cpt++;
    }
    
    val /=double(cpt);
    var /= double(cpt);
    var= (var - val*val);
    return std::make_pair(val, var);
}


std::vector<const char *> options;

// recherche des options
int find( const char *name )
{
    for(int i= 0; i < int(options.size()); i++)
        if(std::string(options[i]) == std::string(name))
            return i;
    return -1;
}

bool find_value( const char *name, const int n )
{
    int index= find(name);
    return (index >= 0 && index + n < int(options.size()));
}

// cherche un flag, renvoie un bool
bool flag( const char *name )
{
    int index= find(name);
    if(index != -1)
        options.erase(options.begin() + index, options.begin() + index +1);
    return (index != -1);
}

const char *value( const char *name, const char *default_value )
{
    int index= find(name);
    if(index < 0 || index + 1 >= int(options.size()))
        return default_value;
    
    const char *str= options[index +1];
    options.erase(options.begin() + index, options.begin() + index +2);
    return str;
}

int value_int( const char *name, const char *default_value )
{
    const char *str= value(name, default_value);
    
    int v= 0;
    if(sscanf(str, "%d", &v) == 1)
        return v;
    
    // erreur parametre invalide.
    printf("invalid int parameter: %s = '%s'\n", name, str);
    
    exit(1);
    return 0;
}

double value_double( const char *name, const char *default_value )
{
    const char *str= value(name, default_value);
    
    double v= 0;
    if(sscanf(str, "%lf", &v) == 1)
        return v;
    
    // erreur parametre invalide.
    printf("invalid double parameter: %s = '%s'\n", name, str);
    
    exit(1);
    return 0;
}


int main( const int argc, const char **argv )
{
    if(argc == 1)
    {
        printf("usage: %s [--pixel x y --error] output.[pfm|hdr|exr] reference.[pfm|hdr|exr] *.[pfm|hdr|exr]\n", argv[0]);
        return 0;
    }
    
    options= std::vector<const char *>(argv, argv + argc);
    
    bool ref_error= false;
    bool export_pixel= false;
    int xpixel= 0, ypixel= 0;
    if(find_value("--pixel", 2))
    {
        int index= find("--pixel");
        if(sscanf(options[index+1], "%d", &xpixel) != 1)
            return 1;
        if(sscanf(options[index+2], "%d", &ypixel) != 1)
            return 1;
        
        options.erase(options.begin() + index, options.begin() + index +3);
        export_pixel= true;
        
        ref_error= flag("--error");
        if(ref_error)
            printf("using mean error images.\n");
    }

    MAXMSE= value_double("--scale", " 0.002");
    int size= value_int("--size", "3");
    
    //~ for(int i= 0; i < int(options.size()); i++)
        //~ printf("%s ", options[i]);
    //~ printf("\n");
    
    //~ printf("--size %d\n", size);
    //~ printf("--scale %lf\n", MAXMSE);
    //~ printf("--pixel %s %d %d\n", export_pixel ? "true" : "false", xpixel, ypixel);
    //~ return 0;
    
    //
    Image ref;
    if(int(options.size()) > 2)
    {
        if(is_exr_image(options[2]))
            ref= read_image_exr(options[2]);
        else if(is_pfm_image(options[2]))
            ref= read_image_pfm(options[2]);
        else if(is_hdr_image(options[2]))
            ref= read_image_hdr(options[2]);
    }
    if(ref.size() == 0)
    {
        printf("[error] loading ref image '%s'...\n", options[2]);
        return 1;
    }
    
    //
    std::vector<Image> images(options.size() -3);
    std::vector<double> spp(options.size() -3);
#pragma omp parallel for     
    for(int i= 3; i < int(options.size()); i++)
    {
        Image image;
        if(is_exr_image(options[i]))
            image= read_image_exr(options[i]);
        else if(is_pfm_image(options[i]))
            image= read_image_pfm(options[i]);
        else if(is_hdr_image(options[i]))
            image= read_image_hdr(options[i]);
        
        //~ images.push_back(image);
        images[i-3]= image;
        
    #if 0
        int n= -1;
        // suppose que le nom de fichier est de la forme 'scene-sampler-nspp.pfm'
        if(sscanf(options[i], "%*[^-]-%*[^-]-%d.pfm", &n) == 1
        || sscanf(options[i], "%*[^-]-%*[^-]-%d.exr", &n) == 1
        || sscanf(options[i], "%*[^-]-%*[^-]-%d.hdr", &n) == 1)
        {
            spp.push_back(std::log10(double(n)));
            printf("    %dspp\n", n);
            //~ printf("%s %dx%d %dspp\n", options[i], images.back().width(),images.back().height(), n);
            //~ printf("%dx%d %dspp ", images.back().width(),images.back().height(), n);
        }
        else
        {
            printf("\n[error] filename '%s' doesn't match scene-sampler-nspp.[pfm|hdr|exr]\n", options[i]);
            exit(1);
        }
    #else
        
        int n= -1;
        // suppose que le nom de fichier est de la forme 'prefix-nspp.[pfm|exr|hdr]'
        const char *ext= strrchr(options[i], '.');
        if(ext)
        {
            if(strcmp(ext, ".pfm") == 0 || strcmp(ext, ".hdr") == 0 || strcmp(ext, ".exr") == 0)
            {
                const char *last= strrchr(options[i], '-');
                if(sscanf(last, "-%d.", &n) == 1)
                {
                    //~ spp.push_back(std::log10(double(n)));
                    spp[i-3]= std::log10(double(n));
                    printf("    %dspp\n", n);
                }
            }
        }
        
        if(n == -1)
        {
            printf("\n[error] filename '%s'... doesn't match 'prefix-nspp.[pfm|hdr|exr]'\n", options[i]);
            exit(1);
        }
    #endif
    }
    printf("done.\n");
    
    Image mse_image(ref.width(), ref.height(), Black());
    Image offset_image(ref.width(), ref.height(), Black());
    
    std::vector<Image> error_images(images.size(), mse_image);
    
#pragma omp parallel for 
    for(int y= size; y < ref.height() - size; y++)
    {
        for(int x= size; x < ref.width() - size; x++)
        {
            std::vector<double> mse;
            std::vector<double> vmin;
            std::vector<double> vmax;
            std::vector<double> var;
            for(int i= 0; i < int(images.size()); i++)
            {
                double v0, v1;
                //~ auto stat= MSElocal(ref, images[i], x - size, y - size, x + size, y + size, v0, v1);        // MSE
                std::pair<double, double> stat;
                
                if(ref_error)
                    stat= Mlocal_error(images[i], x - size, y - size, x + size, y + size, v0, v1);      // Mean Error
                else
                    stat= Mlocal(ref, images[i], x - size, y - size, x + size, y + size, v0, v1);      // Mean
                
                mse.push_back(std::log10(stat.first));
                var.push_back(stat.second);
                vmin.push_back(v0);
                vmax.push_back(v1);
            }
            
            double slope, b;
            llsq(spp.size(), spp.data(), mse.data(), slope, b);
            //~ mse_image(x, y)= colormap(std::abs(slope), 2.0);
            //~ offset_image(x, y)= colormap(std::abs(std::pow(10.0, b)), MAXOFFSET);
            mse_image(x, y)= Color32(std::abs(slope));
            //~ offset_image(x, y)= Color32(std::abs(std::pow(10.0, b)));
            
            {
                // integre l'erreur sur les n images
                std::vector<double> z;
                std::vector<double> f;
                assert(spp.size() == mse.size());
                for(int i= 0; i < int(mse.size()); i++)
                {
                    z.push_back(std::pow(10.0, spp[i]));
                    f.push_back(std::pow(10.0, mse[i]));
                }
                
                int n= 0;
                double v= 0;
                for(int i= 0; i +1 < int(z.size()); i++)
                {
                    double a= z[i];
                    double b= z[i+1];
                    v= v + (b - a) * (f[i] + f[i+1]) / 2;
                    n++;
                }
                v= v / double(n);
                offset_image(x, y)= Color32(v);
            }
            
            if(export_pixel && x == xpixel && y == ypixel)
            {
                char prefix[1024];
                sscanf(options[1], "%[^.].", prefix);
                
                char tmp[1024];
                sprintf(tmp, "pixel_%d_%d_", xpixel, ypixel);
                std::string exportfile= change_prefix(tmp, change_extension(options[1], ".txt"));
                printf("export pixel '%s'...\n\n", exportfile.c_str());
                FILE *out= fopen(exportfile.c_str(), "wt");
                
                fprintf(out, "#Nbpts\t#Mean\t#Fit\t#Min\t#Max\t#Var\n");
                
                //~ fprintf(out, "%0.15lf %0.15lf %0.15lf\n", 0.0, std::pow(10.0, mse[0]), b);
                for(int i= 0; i < int(mse.size()); i++)
                    fprintf(out, "%0.15lf\t%0.15lf\t%0.15lf\t%0.15lf\t%0.15lf\t%0.15lf\n", std::pow(10.0, spp[i]), std::pow(10.0, mse[i]), std::pow(10.0, slope * spp[i] + b), vmin[i], vmax[i], var[i]);
                    //            spp     mse     fit     vmin    vmax    var
                
                fclose(out);
            }
            
            for(int i= 0; i < int(images.size()); i++)
            {
                Color32 d= ref(x, y) - images[i](x, y);
                double v= (d.r*d.r + d.g*d.g + d.b*d.b) / 3;
                
                //~ Color32 color= colormap(v, MAXMSE / std::pow(10.0, spp[i]));
                Color32 color= Color32(v);
                error_images[i](x, y)= color;
            }
        }
    }
    
    if(export_pixel)
    {
        for(int y= ypixel - size; y <= ypixel + size; y++)
        for(int x= xpixel - size; x <= xpixel + size; x++)
            mse_image(x, y)= Black();

        char tmp[1024];
        sprintf(tmp, "pixel_%d_%d_", xpixel, ypixel);
        write_image_hdr(mse_image, ::change_prefix(tmp, ::change_extension(options[1], ".hdr")).c_str());
        return 0;
    }
    
    if(is_exr_image(options[1]))
        write_image_exr(mse_image, options[1]);
    else if(is_pfm_image(options[1]))
        write_image_pfm(mse_image, options[1]);
    else if(is_hdr_image(options[1]))
        write_image_hdr(mse_image, options[1]);
    
    // ajoute une echelle a droite de l'image
    //~ write_image_exr(mse_image, ::change_prefix("slopes_", ::change_extension(options[1], ".exr")).c_str());
    for(int y= 0; y < mse_image.height(); y++)
    {
        for(int x= 0; x < mse_image.width(); x++)
        {
            float v= mse_image(x, y).r;
            mse_image(x, y)= colormap(v, 2.0);
        }
        
        //~ for(int x=  mse_image.width() -32; x <  mse_image.width() -16; x++)
            //~ mse_image(x, y)= Color32(2.0 * double(y) / double(mse_image.height() -1));
        
        for(int x=  mse_image.width() -16; x <  mse_image.width(); x++)
            mse_image(x, y)= colormap(2.0 * double(y) / double(mse_image.height() -1), 2.0);
    }
    //~ write_image(mse_image, change_prefix("slopes_", change_extension(options[1], ".png")).c_str());
    
    //~ // ajoute une echelle a droite de l'image
    //~ write_image_exr(offset_image, change_prefix("offset_", change_extension(options[1], ".exr")).c_str());
    
    //~ for(int y= 0; y < offset_image.height(); y++)
    //~ {
        //~ for(int x= 0; x < offset_image.width(); x++)
        //~ {
            //~ float v= offset_image(x, y).r;
            //~ offset_image(x, y)= colormap(v, MAXOFFSET);
        //~ }
        
        //~ for(int x=  offset_image.width() -32; x <  offset_image.width() -16; x++)
            //~ offset_image(x, y)= Color32(MAXOFFSET * double(y) / double(offset_image.height() -1));
        
        //~ for(int x=  mse_image.width() -16; x <  mse_image.width(); x++)
            //~ offset_image(x, y)= colormap(MAXOFFSET * double(y) / double(offset_image.height() -1), MAXOFFSET);
    //~ }
    //~ write_image(offset_image, change_prefix("offset_", change_extension(options[1], ".png")).c_str());
    
    // ajoute une echelle a droite de l'image
    printf("using scale %lf...\n", MAXMSE);
#pragma omp parallel for
    for(int i= 0; i < int(images.size()); i++)
    {
        //~ write_image_exr(error_images[i], change_prefix("error_", change_extension(options[i+3], ".exr")).c_str());
        
        for(int y= 0; y < error_images[i].height(); y++)
        {
            for(int x= 0; x < error_images[i].width(); x++)
            {
                float v= error_images[i](x, y).r;
                error_images[i](x, y)= colormap(v, MAXMSE / std::pow(10.0, spp[i]));
            }
        }
        
        for(int y= 0; y < error_images[i].height(); y++)
        {
            //~ for(int x= error_images[i].width() -32; x < error_images[i].width() -16; x++)
                //~ error_images[i](x, y)= Color32(MAXMSE / std::pow(10.0, spp[i]) * double(y) / double(error_images[i].height() -1));
            
            for(int x= error_images[i].width() -16; x < error_images[i].width(); x++)
                error_images[i](x, y)= colormap(MAXMSE / std::pow(10.0, spp[i]) * double(y) / double(error_images[i].height() -1), MAXMSE / std::pow(10.0, spp[i]));
        }
        write_image(error_images[i], change_prefix("error_", change_extension(options[i+3], ".png")).c_str());
    }
    
    printf("\n");
    
    return 0;
}
