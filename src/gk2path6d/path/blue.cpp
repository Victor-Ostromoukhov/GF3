
#include <cstdarg>
#include <cassert>
#include <climits>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <algorithm>
#include <omp.h>

#include "options.h"

#include "real.h"
#include "vec.h"
#include "mat.h"
#include "color.h"

#include "mesh.h"
#include "wavefront.h"
#include "mesh_cache.h"

#include "image.h"
#include "image_hdr.h"
#include "image_exr.h"

#include "ray.h"
#include "triangle.h"

#include "real.h"

#include "framebuffer.h"

#include "sampler.h"
#include "sampler_sobol.h"
#include "sampler_sobolpp.h"
#include "sampler_nd.h"
#include "sampler_s19.h"

#include "scene.h"

#include "points.h"
#include "brdf.h"
#include "orbiter.h"
#include "camera.h"


// utilitaire
struct Format
{
    char text[1024];
    
    Format( const char *_format, ... )
    {
        text[0]= 0;     // chaine vide
        
        // recupere la liste d'arguments supplementaires
        va_list args;
        va_start(args, _format);
        vsnprintf(text, sizeof(text), _format, args);
        va_end(args);
    }
    
    ~Format( ) {}
    
    // conversion implicite de l'objet en chaine de caracteres stantard
    operator const char *( ) { return text; }
    
    operator std::string( ) { return std::string(text); }
};


/*! renvoie le nom d'un fichier.
    basename("path/to/file") == "file"
    basename("file") == "file"
 */
static
std::string file_basename( const char *filename )
{
    std::string path= filename;
#ifndef WIN32
    std::replace(path.begin(), path.end(), '\\', '/');   // linux, macos : remplace les \ par /.
    size_t slash = path.find_last_of( '/' );
    if(slash != std::string::npos)
        return path.substr(slash +1, std::string::npos);
    else
        return filename;
#else
    std::replace(path.begin(), path.end(), '/', '\\');   // windows : remplace les / par \.
    size_t slash = path.find_last_of( '\\' );
    if(slash != std::string::npos)
        return path.substr(slash +1, std::string::npos);
    else
        return filename;
#endif
}


// renderer data
std::vector<Source> sources;

Scene scene;
Camera camera;
Options options;


// renderer techniques

// simple direct lighting
Color direct( Scene& scene, const Float t, 
    const UniformSource& points, Sampler& sampler, 
    const Raydata& primary )
{
    Point o= primary.o();
    Point p= primary.hitp();
    Vector pn= scene.normal(primary);    
    Vector po=normalize(Vector(p, o));
    
    Color sample;
    
    // choisit un point sur une source
    Float u1= 0;
    Float u2= 0;
    Float u3= 0;
    if(options.remap_bounce3d)
    {
        // utilise les memes dimensions pour sampler les 2 strategies / 2 integrales
        u1= sampler.sample1();
        u2= sampler.sample2();
        u3= sampler.sample3();
        if(options.remap1)
            u1= 0;  // force le selecteur a 0, choisit toujours le premier terme
    }
    else if(options.remap_bounce2d)
    {
        // utilise les memes dimensions pour sampler les 2 strategies / 2 integrales
        u1= sampler.sample1();
        u2= sampler.sample2();
    }
    
    Float pdf;
    LightSample light;
    if(options.remap_bounce3d)
        light= points.sample(u1, u2, u3, pdf);
    else if(options.remap_bounce2d)
        light= points.sample(u1, u2, pdf);
    else
        light= points.sample(sampler, pdf);
    
    Vector l= normalize(Vector(p, light.s));
    Float cos_theta= std::max(Float(0), dot(pn, l));
    Float cos_theta_s= std::max(Float(0), -dot(light.sn, l));
    Float g= cos_theta_s / distance2(p, light.s);
    
    if(cos_theta * g > 0 
    && scene.visible(t, p + pn * Float(.001), light.s + light.sn * Float(.001)))
    {
        int lod= 0;
        Color diffuse= scene.diffuse_texture(primary, lod);
        const Material& material= scene.material(primary);
        
        Float d= material.diffuse.power();
        Float s= material.specular.power();
        if(options.force_diffuse_material) { d= 1; s= 0; }      // force matiere diffuse
        Float kd= d / (d+s);
        Float ks= s / (d+s);
        
        Blinn distribution(material.ns);
        Brdf brdf(kd, diffuse * material.diffuse, ks, distribution, material.specular);
        
        World world(pn);
        Color fr= brdf.f(world.inverse(l), world.inverse(po));
        
        sample= sample + light.emission * fr * cos_theta * g / pdf;
    }
    
    return sample;
}


// direct lighting : multiple importance sampling 
Color direct( Scene& scene, const Float t, 
    const UniformSource& points, Sampler& sampler, 
    const Point& o, const Point& p, const Vector& pn, Brdf& brdf, Raydata& indirect_hit )
{
    World world(pn);
    Vector po=normalize(Vector(p, o));
    
    Color sample;
    Float u1= 0;
    Float u2= 0;
    Float u3= 0;
    if(options.remap_bounce3d)
    {
        // utilise les memes dimensions pour sampler les 2 strategies / 2 integrales
        u1= sampler.sample1();
        u2= sampler.sample2();
        u3= sampler.sample3();
        if(options.remap1)
            u1= 0;  // force le selecteur a 0, choisit toujours le premier terme
    }
    else if(options.remap_bounce2d)
    {
        // utilise les memes dimensions pour sampler les 2 strategies / 2 integrales
        u1= sampler.sample1();
        u2= sampler.sample2();
    }
    
    // mis : strategie 1
    // choisit un point sur une source
    Float light_pdf;
    {
        LightSample light;
        if(options.remap_bounce3d)
            light= points.sample(u1, u2, u3, light_pdf);
        else if(options.remap_bounce2d)
            light= points.sample(u1, u2, light_pdf);
        else
            light= points.sample(sampler, light_pdf);
        assert(light_pdf > 0);
        
        Vector l= normalize(Vector(p, light.s));
        Float cos_theta= std::max(Float(0), dot(pn, l));
        Float cos_theta_s= std::max(Float(0), -dot(light.sn, l));
        Float g= cos_theta_s / distance2(p, light.s);
        Color fr= brdf.f(world.inverse(l), world.inverse(po));
        
        if(cos_theta * g > 0 
        && scene.visible(t, p + pn * Float(.001), light.s + light.sn * Float(.001) ))
        {
            // evalue l'autre strategie
            Float pdf2= brdf.eval(world.inverse(l), world.inverse(po)) * g;
            
            // poids
            //~ Float w= light_pdf / (light_pdf + pdf2);        // balance heuristic
            Float w= (light_pdf*light_pdf) / (light_pdf*light_pdf + pdf2*pdf2);        // power heuristic
            sample= sample + w * light.emission * fr * cos_theta * g / light_pdf;
        }
    }

    // mis : strategie 2
    // choisit une direction 
    Float brdf_pdf;
    {
        Vector l;
        if(options.remap_bounce3d)
            l= world(brdf.sample(u1, u2, u3, world.inverse(po), brdf_pdf));
        else if(options.remap_bounce2d)
            l= world(brdf.sample(u1, u2, world.inverse(po), brdf_pdf));
        else
            l= world(brdf.sample(sampler, world.inverse(po), brdf_pdf));
        
        if(brdf_pdf > 0)
        {
            Raydata shadow(p + pn * Float(.001), l, t);
            if(scene.intersect(shadow))
            {
                Vector shadow_hitn= scene.normal(shadow);
                
                const Material& shadow_material= scene.material(shadow);
                if(dot(shadow_hitn, l) * shadow_material.emission.power() < 0)
                {
                    Float cos_theta= std::max(Float(0), dot(pn, l));
                    Float cos_theta_s= std::max(Float(0), -dot(shadow_hitn, l));
                    Float g= cos_theta_s / distance2(p, shadow.hitp());
                    Color fr= brdf.f(world.inverse(l), world.inverse(po));

                    // evalue l'autre strategie
                    brdf_pdf= brdf_pdf * g;
                    Float pdf2= points.eval(shadow.hitp(), shadow_material);
                    
                    // poids
                    //~ Float w= brdf_pdf / (brdf_pdf + pdf2);    // balance heuristic
                    Float w= (brdf_pdf*brdf_pdf) / (brdf_pdf*brdf_pdf + pdf2*pdf2);    // power heuristic
                    sample= sample + w * shadow_material.emission * fr * cos_theta * g / brdf_pdf;
                }
                else if(shadow_material.emission.power() == 0)
                    // re-utilise l'echantillon pour prolonger le chemin
                    // sauf sur une source de lumiere...
                    indirect_hit= shadow;
            }
        }
    }
    
#if 0
    // version directe multiple importance sampling, regenere une direction + rayon
    // cf pbrt...
    #ifndef REMAP_BOUNCE
    {
        Float pdf;
    #if defined REMAP_BOUNCE_3D
        Vector l= world(brdf.sample(u1, u2, u3, world.inverse(po), pdf));
    #elif defined REMAP_BOUNCE_2D
        Vector l= world(brdf.sample(u1, u2, world.inverse(po), pdf));
    #else
        Vector l= world(brdf.sample(sampler, world.inverse(po), pdf));
    #endif
        
        if(pdf > 0)
        {
            Raydata shadow(p + pn * Float(.001), l, t);
            if(scene.intersect(shadow))
            {
                const Material& shadow_material= scene.material(shadow);
                if(shadow_material.emission.power() == 0)
                    // re-utilise l'echantillon pour prolonger le chemin
                    // sauf sur une source de lumiere...
                    indirect_hit= shadow;
            }
        }
    }
    #endif
#endif

    return sample;
}


Color path( Scene& scene, const Float t,
    const int depth, const UniformSource& points, Sampler& sampler, 
    const Raydata& primary_hit )
{
    Color sample= Black();

    if(!primary_hit.intersect())
        return sample;
    
    Point o= primary_hit.o();
    Raydata hit= primary_hit;
    
    Color weight= White();
    for(int i= 0; i <= depth; i++)
    {
        int lod= (i == 0) ? 0 : 16;
        Color diffuse= scene.diffuse_texture(hit, lod);
        const Material& material= scene.material(hit);
        assert(material.emission.power() == 0);
        
        Float d= material.diffuse.power();
        Float s= material.specular.power();
        if(options.force_diffuse_material) { d= 1; s= 0; }      // force matiere diffuse
        Float kd= d / (d+s);
        Float ks= s / (d+s);
        
        Blinn distribution(material.ns);
        Brdf brdf(kd, diffuse * material.diffuse, ks, distribution, material.specular);
        
        // estime le direct en p, terme L_(i)
        Point p= hit.hitp();
        Vector po= normalize(Vector(p, o));
        Vector n= scene.normal(hit);
        if(dot(po, n) < 0)
            n= -n;
        
        Raydata indirect_hit;
        Color direct_sample= direct(scene, t, points, sampler, o, p, n, brdf, indirect_hit);
        
        sample= sample + weight * direct_sample;
        
        // prepare l'estimation de l'indirect en p, terme L_(i+1)
        if(!indirect_hit.intersect())
            break;
        
        Point q= indirect_hit.hitp();
        Vector l= normalize(Vector(p, q));
        
        World world(n);
        Float cos_theta= std::max(Float(0), dot(n, l));
        Color fr= brdf.f(world.inverse(l), world.inverse(po)) ;
        Float pdf= brdf.eval(world.inverse(l), world.inverse(po));
        
        weight= weight * fr * cos_theta / pdf;
    #if 0
	// roulette russe, version pbrt3
        if(i > 3)
        {
            Float q= std::max(Float(.05), 1 - weight.power());
            if(sampler.sample1() < q)
                return Color(sample, i);
            
            weight= weight / (1 - q);
        }
    #endif
        
        o= p;
        hit= indirect_hit;
    }
    
    return sample;
}


Color ambient( Scene& scene, const Float t, 
    Sampler& sampler, Direction *direction, const Raydata& primary )
{
    assert(direction != nullptr);
    
    Color sample= Black();
    
    Point o= primary.o();
    Point p= primary.hitp();
    Vector pn= scene.normal(primary);
    Vector po= normalize(Vector(p, o));
    if(dot(po, pn) < 0)
        pn= -pn;
    
#if 1
    Vector wl= direction->sample(sampler.sample1(), sampler.sample2());
    Float pdf= direction->eval(wl);
    
#else
    // genere une direction sur une spirale de fibonacci
    Fibonacci dd(sampler.size(), sampler.seed());
    Vector wl= dd.direction(sampler.index());
    Float pdf= 1 / Float(M_PI);
#endif
    
    // ciel visible ?
    World world(pn);
    Vector w= world(wl);
    Raydata shadow(p + Float(.001) * pn, w, t);
    
    if(!scene.intersect(shadow))
    {
        int lod= 0;
        Color diffuse= scene.diffuse_texture(primary, lod);
        
        Float cos_theta= std::max(Float(0), dot(pn, w));
        //~ sample= sample + Color(Float(1) / Float(M_PI)) * cos_theta / pdf;
        sample= sample + diffuse / Float(M_PI) * cos_theta / pdf;
    }
    
    return sample;
}


// predicat pour trier des couleurs
bool color_grey_less( const Color& a, const Color& b )
{
    //~ return (a.r + a.g + a.b) < (b.r +b.g + b.b);
    
    // https://en.wikipedia.org/wiki/Grayscale
    Color y= Color(Float(0.299), Float(0.587), Float(0.114));           // rec601
    //~ Color y= Color(Float(0.2126), Float(7152), Float(0.0722));          // ITU-R BT.709 standard used for HDTV television
    //~ Color y= Color(Float(0.2627), Float(0.6780), Float(0.0593));        // ITU-R BT.2100 standard for HDR television
    return (y.r * a.r + y.g * a.g + y.b * a.b) < (y.r * b.r + y.g * b.g + y.b * b.b);
}

struct color_index_less
{
    const std::vector<Color>& values;
    color_index_less( const std::vector<Color>& v ) : values(v) {}
    
    bool operator() ( const int a, const int b ) const
    {
        return color_grey_less(values[a], values[b]);
    }
};

int main( int argc, char *argv[] )
{
    options= Options((const char **) argv, (const char **) argv + argc);
    
    // charge le mesh et extrait toutes les donnees
    {
        if(options.scene_filename)
            // charge une description de scene, si possible
            scene= read_scene(options.scene_filename);
        else if(options.file_filename)
            // sinon, fabrique une scene avec un seul objet statique
            scene.attach( SceneObject(options.file_filename) );
        else
        {
            printf("[error] no scene...\n");
            return 1;
        }
        
        scene.build();  // construit le stbvh
        
        // ensemble de sources de lumieres
        scene.build_sources(sources);
        
        if(sources.empty() && !options.force_ambient_occlusion)
        {
            printf("[error] no light sources...\n");
            return 1;
        }
    }
    
    // charge le point de vue
    Orbiter orbiter;
    {
        const char *filename= nullptr;
        if(options.orbiter_filename)
            filename= options.orbiter_filename;
        else 
            filename= scene.camera_filename.c_str();
        
        if(filename == nullptr || orbiter.read_orbiter(filename) < 0)
        {
            printf("[error] no camera...\n");
            return 1;
        }
    }
    
    Lens lens;
    if(options.sample_lens)
    {
        const char *filename= nullptr;
        if(options.lens_filename)
            filename= options.lens_filename;
        else 
            filename= scene.lens_filename.c_str();
        
        if(filename == nullptr || lens.read_lens(filename) < 0)
        {
            printf("[error] no lens...\n");
            return 1;
        }
    }
    
    // 
    camera= Camera(options.width, options.height, 45, orbiter, lens);
    Framebuffer image(options.width, options.height);

    // distribution de sources
    UniformSource points(sources);
    printf("%d sources.\n", points.size());
    
    // masque blue noise
    Image mask= read_image_pfm("data/mask.pfm");
    if(mask.size() == 0)
    {
        printf("[error] no blue noise mask...\n");
        return 1;
    }
    
    // samplers 
    Sampler *main_sampler= nullptr;
    
    size_t main_samples= size_t(options.width) * size_t(options.height) * size_t(options.samples_per_pixel);
    size_t pixel_samples= options.samples_per_pixel;
    int width= options.width;
    int height= options.height;
    
    if(main_samples == 0)
    {
        main_samples= options.samples;
        pixel_samples= options.samples / (width * height);
    }
    
    assert(main_samples != 0);
    printf("using %ld samples: ", main_samples);
    printf("%dx%d %dspp\n", width, height, int(pixel_samples));
    
    //~ unsigned int seed= seed_generator();
    unsigned int seed= 0;
    //~ main_sampler= new OwenPixelSampler(options.dimensions_remap, pixel_samples, options.sobol_filename, options.owen_bits, /* seed */ seed);
    main_sampler= new CSobolPixelSampler(options.dimensions, pixel_samples, options.sobol_filename, options.owen_bits, /* seed */ seed);
    assert(main_sampler != nullptr);
    
    if(options.samples_remap_string)
        main_sampler->remap_dimensions(options.samples_remap_dimensions, options.samples_remap);
    else
        main_sampler->remap_dimensions(options.dimensions, options.dimensions_remap);
    
    // generateur de directions pour l'ambient
    Direction *ambient_direction= nullptr;
    if(options.force_ambient_occlusion || points.size() == 0)
    {
        if(options.ambient_use_disk_uniform)
            ambient_direction= new ConcentricDiskUniformDirection();
        else if(options.ambient_use_disk_cosine)
            ambient_direction= new ConcentricDiskCosineDirection();
        else if(options.ambient_use_uniform)
            ambient_direction= new UniformDirection();
        else
            ambient_direction= new CosineDirection();
    }
    
    //
    std::chrono::high_resolution_clock::time_point cpu_start= std::chrono::high_resolution_clock::now();
    
//~ #undef THREADS
#define THREADS
    
    int nthreads= 1;
#ifdef THREADS // debug, pas de parallelisation
    if(options.threads == 0)
        nthreads= omp_get_max_threads();
    else
        nthreads= std::min(options.threads, omp_get_max_threads());
    
    if(options.debug_pixel_filename && options.debug_export_samples)
        nthreads= 1;    // force un seul thread...
        
    omp_set_num_threads(nthreads);
#endif
    
    printf("%d threads.\n", nthreads);
    
    // initialise un sampler par thread
    std::vector<Sampler *> samplers;
    for(int i= 0; i < nthreads; i++)
        samplers.push_back(main_sampler->clone(seed));
    
    // seeds par pixel
    BlueTile image_seeds(options.width, options.dimensions, 1);
    assert(options.width == options.height);
    
    int lines= 0;
#pragma omp parallel for schedule (dynamic, 1)
    for(int py=0; py < image.height(); py++)
    {
    #pragma omp atomic
        lines++;
        
        if(100*lines / image.height() != 100*(lines+1) / image.height())
        {
            printf("%d/%d %d%%...\r", lines, image.height(), 100*lines / image.height());
            fflush(stdout);
        }
        
        // recupere le sampler associe au thread
        unsigned int tid= omp_get_thread_num();
        assert(tid < samplers.size());
        Sampler& sampler= *samplers[tid];
        
        std::random_device seed_generator;
        
        constexpr int N= 256;
        std::vector<Color> values(N);
        std::vector<int> index(N);
        std::vector<uint32_t> seeds(options.dimensions*N);
        
        for(int px= 0; px < image.width(); px++)
        {
            for(int v= 0; v < values.size(); v++)
            {
            #if 0
                unsigned int seed= seed_generator();
                sampler.seed(seed);
            #else
                for(int d= 0; d < options.dimensions; d++)
                    seeds[v*options.dimensions + d]= seed_generator();
                
                sampler.pixel_seeds(options.dimensions, seeds.data() + v*options.dimensions);
            #endif
                
                Color value;
                for(int s= 0; s < pixel_samples; s++)
                {
                    sampler.index(s);
                    
                    Float x= Float(px) + Float(.5);
                    Float y= Float(py) + Float(.5);
                    if(options.sample_subpixel)
                    {
                        sampler.dimension(0);       // force les dimensions utilisees pour generer le rayon primaire
                        x= Float(px) + sampler.sample1();
                        y= Float(py) + sampler.sample2();
                    }
                    
                    Float lens_u= Float(0.5);
                    Float lens_v= Float(0.5);
                    if(options.sample_lens) 
                    { 
                        sampler.dimension(2);    // force les dimensions utilisees pour generer le rayon primaire
                        lens_u= sampler.sample3(); 
                        lens_v= sampler.sample4();
                    }
                    
                    Float time= 0;
                    if(options.sample_time)
                    {
                        sampler.dimension(4);    // force les dimensions utilisees pour generer le rayon primaire
                        time= sampler.sample5();
                    }
                    
                    // force les dimensions utilisees pour l'eclairage direct et les rebonds
                    sampler.dimension(5);
                    
                    Point o;
                    Vector d;
                    Float lens_pdf= camera.sample(x, y, lens_u, lens_v, o, d);
                    
                    Raydata ray(o, d, time);
                    Color sample;
                    if(scene.intersect(ray))
                    {
                        Vector hitn= scene.normal(ray); // normale interpolee du triangle
                        const Material& material= scene.material(ray);
                        
                        // accumule l'emission de la matiere, suppose que les normales sont correctes... 
                        // evite de creer des sources qui emettent sur les 2 faces...
                        if(dot(hitn, ray.d()) < 0)
                            sample= sample + material.emission / lens_pdf;
                        
                        if(material.emission.power() == 0)
                        {
                            if(!options.force_ambient_occlusion && points.size() > 0)
                            {
                                if(options.force_direct_lighting)
                                    // direct lighting
                                    sample= sample + direct(scene, time, points, sampler, ray) / lens_pdf;
                                else
                                    // path tracing
                                    sample= sample + path(scene, time, /*depth */ options.path_length, points, sampler, ray) / lens_pdf;
                            }
                            else
                            {
                                // pas de sources, ambient occlusion
                                sample= sample + ambient(scene, time, sampler, ambient_direction, ray) / lens_pdf;
                            }
                        }
                    }
                    
                    value= value + sample;
                }
                
                values[v]= value / Float(pixel_samples);
            }
            
            // trier les valeurs
            for(int v= 0; v < values.size(); v++)
                index[v]= v;
            //~ std::sort(values.begin(), values.end(), color_grey_less);
            std::sort(index.begin(), index.end(), color_index_less(values));
            
            Color32 m= mask(px % mask.width(), py % mask.height());
            Float g= (Float(m.r) + Float(m.g) + Float(m.b)) / 3;        // cf mask.pfm
            //~ Color value= values[g*int(values.size())];
            
            int seed_id= index[g*int(values.size())];
            
            Color value= values[seed_id];
            image.splat(px, py, Color(value, 1));
            
            // conserve les seeds
            for(int d= 0; d < options.dimensions; d++)
                image_seeds(px, py, d, 0)= seeds[seed_id*options.dimensions + d];
        }
    }
    printf("\ndone.\n");
    
    // prepare les samplers pour le prochain run
    for(int tid= 0; tid < int(samplers.size()); tid++)
        delete samplers[tid];

    // affiche quelques stats
    std::chrono::high_resolution_clock::time_point cpu_stop= std::chrono::high_resolution_clock::now();
    int cpu_time= std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
    printf("cpu  %ds %03dms\n", (int) (cpu_time / 1000), (int) (cpu_time % 1000));
    printf("%fMsamples/s\n", (image.width() * image.height() * options.samples_per_pixel / Float(1000000.0)) / (Float(cpu_time) / Float(1000.0)));
    
    // enregistre les images
    write_image_pfm(image.flatten(), Format("%s-%05lu.pfm", options.output_prefix, pixel_samples));
//    printf("writing output image '%s'...\n", options.output_filename);
//    write_image(image.flatten(2.2), options.output_filename);
    
    image_seeds.write(Format("%s_seeds_%dd-%05lu.dat", options.output_prefix, options.dimensions, pixel_samples));
    
    {
        SeedTile tmp(image_seeds.width(), image_seeds.dimensions());
        for(int y= 0; y < image_seeds.height(); y++)
        for(int x= 0; x < image_seeds.width(); x++)
        for(int d= 0; d < image_seeds.dimensions(); d++)
            tmp(x, y, d)= image_seeds(x, y, d, 0);
        tmp.saveTile("seedmap.dat");
    }
    
    main_sampler->release();
    delete main_sampler;

    return 0;
}
