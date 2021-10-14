
#include <cstdarg>
#include <cassert>
#include <climits>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <algorithm>
#include <omp.h>

#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <eigen3/Eigen/Dense>

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
#include "sampler_zsampler.h"
#include "sampler_morton.h"

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

std::vector<size_t> path_depth; // total samples
std::vector<int> last_depth;    // pixel samples

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
    {
        indirect_hit= Raydata();
        
        Float pdf;
        Vector l;
        if(options.remap_bounce3d)
            l= world(brdf.sample(u1, u2, u3, world.inverse(po), pdf));
        else if(options.remap_bounce2d)
            l= world(brdf.sample(u1, u2, world.inverse(po), pdf));
        else
            l= world(brdf.sample(sampler, world.inverse(po), pdf));
            
        if(pdf > 0)
        {
            Raydata shadow(p + pn * Float(.001), l, t);
            if(scene.intersect(shadow))
            {
                Vector shadow_hitn= scene.normal(shadow);
                const Material& shadow_material= scene.material(shadow);
                if(shadow_material.emission.power() == 0)
                    indirect_hit= shadow;
            }
        }
    }
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
    
    //~ std::vector<int> dimensions;
    
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
        
        //~ {
            //~ static int done= 0;
            
            //~ if(!done)
            //~ {
                //~ dimensions.push_back( sampler.dimension() );
                
                //~ if(i == depth)
                //~ {
                    //~ done= 1;
                    
                    //~ // extraction des dimensions utilisees pour le premier sample...
                    //~ for(int d= 0; d < int(dimensions.size()); d++)
                        //~ printf("index %lu depth %d: dimension %d\n", sampler.index(), d, dimensions[d]);
                //~ }
            //~ }
        //~ }
        
        Raydata indirect_hit;
        Color direct_sample= direct(scene, t, points, sampler, o, p, n, brdf, indirect_hit);

        if(i == depth)
            // n'accumule le direct que sur le dernier rebond indirect
            sample= sample + weight * direct_sample;
        
        //~ // comptabilise les chemins existants par dimension
        //~ if(direct_sample.grey() > 0)
        //~ {
            //~ assert(i < int(path_depth.size()));
            //~ #pragma omp atomic
            //~ path_depth[i]++;
        //~ }
        
        //~ // export 4d functions
        //~ if(i == 0 && sampler.dimensions() == 4)
            //~ sample= sample + weight * direct_sample;
        //~ else if(i == 1 && sampler.dimensions() == 6)
            //~ sample= sample + weight * direct_sample;
        
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



Color path1_rank( Scene& scene, const Float t,
    const UniformSource& points, Sampler& sampler, 
    const Raydata& primary )
{
    Color sample= Black();
    if(!primary.intersect())
        return sample;
    
    // get samples 
    const int rankN= options.samples_per_pixel;
    
    std::vector<Float> samples1(rankN);
    std::vector<Float> samples2(rankN);
    std::vector<Float> samples3(rankN);
    std::vector<Float> samples4(rankN);
    
#if 0
    // utilise le sampler global
    size_t sampler_index= sampler.index();
    for(int i= 0; i < rankN; i++)
    {
        sampler.index(sampler_index + i);
        sampler.dimension(5);
        
        samples1[i]= sampler.sample1();
        samples2[i]= sampler.sample2();
        samples3[i]= sampler.sample3();
        samples4[i]= sampler.sample4();
    }
#else
    // genere une grille 2D+2D
    // + rotation pour chaque pixel
    Float rx= sampler.sample1();
    Float ry= sampler.sample2();
    
    int n= std::sqrt(rankN);
    assert(n*n == rankN);
    for(int i= 0; i < rankN; i++)
    {
        Float x= (i % n) / Float(n) + rx / Float(n);
        Float y= (i / n) / Float(n) + ry / Float(n);
        
        samples1[i]= x;
        samples2[i]= y;
        samples3[i]= x;
        samples4[i]= y;
    }
#endif

    // prepare NxN matrix 
    Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> matrix(rankN, rankN);
    matrix.setZero();
    
    for(int i= 0; i < rankN; i++)
    {
        // brdf primaire
        int lod= 16;
        Color diffuse= scene.diffuse_texture(primary, lod);
        const Material& material= scene.material(primary);
        assert(material.emission.power() == 0);
        
        Float d= material.diffuse.power();
        Float s= material.specular.power();
        if(options.force_diffuse_material) { d= 1; s= 0; }      // force matiere diffuse
        Float kd= d / (d+s);
        Float ks= s / (d+s);
        
        Blinn distribution(material.ns);
        Brdf brdf(kd, diffuse * material.diffuse, ks, distribution, material.specular);
        
        Point o= primary.o();
        Point p= primary.hitp();
        Vector po= normalize(Vector(p, o));
        Vector pn= scene.normal(primary);
        if(dot(po, pn) < 0)
            pn= -pn;
            
        World world(pn);
        
        // rebond
        //~ Float u1= sampler.sample1();
        //~ Float u2= sampler.sample2();
        Float u1= samples1[i];
        Float u2= samples2[i];
        
        Float brdf_pdf= 0;
        Vector w= world(brdf.sample(u1, u2, world.inverse(po), brdf_pdf));
        if(brdf_pdf > 0)
        {
            Raydata indirect(p + pn * Float(.001), w, t);
            if(scene.intersect(indirect))
            {
                const Material& material= scene.material(indirect);
                if(material.emission.power() == 0)
                {
                    Point q= indirect.hitp();
                    Vector qp= normalize(Vector(q, p));
                    Vector qn= scene.normal(indirect);
                    if(dot(qp, qn) < 0)
                        qn= -qn;
                        
                    Float cos_theta= std::max(Float(0), dot(pn, w));
                    Color fr= brdf.f(world.inverse(w), world.inverse(po)) ;
                    Float pdf= brdf.eval(world.inverse(w), world.inverse(po));
                    
                    //~ Color weight= fr * cos_theta / pdf;
                    Color weight= fr * cos_theta;
                    
                    // direct
                    for(int j= 0; j < rankN; j++)
                    {
                        //~ Float u3= sampler.sample3();
                        //~ Float u4= sampler.sample4();
                        Float u3= samples3[j];
                        Float u4= samples4[j];
                        
                        Float light_pdf;
                        LightSample light= points.sample(u3, u4, light_pdf);
                        
                        Vector shadow= normalize(Vector(q, light.s));
                        Float cos_theta= std::max(Float(0), dot(qn, shadow));
                        Float cos_theta_s= std::max(Float(0), -dot(light.sn, shadow));
                        Float g= cos_theta_s / distance2(q, light.s);
                        
                        if(cos_theta * g > 0 
                        && scene.visible(t, q + qn * Float(.001), light.s + light.sn * Float(.001)))
                        {
                            int lod= 16;
                            Color diffuse= scene.diffuse_texture(indirect, lod);
                            const Material& material= scene.material(indirect);
                            assert(material.emission.power() == 0);
                            
                            Float d= material.diffuse.power();
                            Float s= material.specular.power();
                            if(options.force_diffuse_material) { d= 1; s= 0; }      // force matiere diffuse
                            Float kd= d / (d+s);
                            Float ks= s / (d+s);
                            
                            Blinn distribution(material.ns);
                            Brdf brdf(kd, diffuse * material.diffuse, ks, distribution, material.specular);
                            
                            World world(qn);
                            Color fr= brdf.f(world.inverse(shadow), world.inverse(qp));
                            
                            //~ sample= sample + light.emission * weight * fr * cos_theta * g / light_pdf;
                            //~ matrix(i, j)= (light.emission * weight * fr * cos_theta * g / light_pdf).power();
                            matrix(i, j)= (light.emission * weight * fr * cos_theta * g).power();
                            //~ sample= sample + (light.emission * weight * fr * cos_theta * g / light_pdf) / Float(rankN*rankN);
                        }
                    }
                }
            }
        }
    }
    
    //~ return sample;
    
    //~ // compute rank
    //~ Eigen::FullPivLU< Eigen::Matrix<Float, rankN, rankN> > LU(matrix);
    //~ LU.setThreshold(1e-5);
    //~ int rank= LU.rank();
    //~ return Color(rank);
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>> solver;
    solver.compute(matrix * matrix.transpose());
    double val= matrix.squaredNorm() -  solver.eigenvalues().reverse()[0];    
    return Color(val);
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
    
    if(options.debug_pixel_filename)
    {
        if(options.debug_pixel_x < 0 || options.debug_pixel_x >= options.width
        || options.debug_pixel_y < 0 || options.debug_pixel_y >= options.height)
        {
            printf("[error] debug pixel: invalid coordinates %d %d\n", options.debug_pixel_x, options.debug_pixel_y);
            return 1;
        }
    }
    
    Image debug_reference;
    if(options.debug_reference)
    {
        if(is_exr_image(options.debug_reference))
            debug_reference= read_image_exr(options.debug_reference);
        else if(is_pfm_image(options.debug_reference))
            debug_reference= read_image_pfm(options.debug_reference);
        else if(is_hdr_image(options.debug_reference))
            debug_reference= read_image_hdr(options.debug_reference);
        
        if(debug_reference.width() != options.width || debug_reference.height() != options.height)
        {
            printf("[error] debug reference: invalid image size %dx%d\n", debug_reference.width(), debug_reference.height());
            return 1;
        }
    }
    
    // 
    camera= Camera(options.width, options.height, 45, orbiter, lens);
    Framebuffer image(options.width, options.height);

    // distribution de sources
    UniformSource points(sources);
    printf("%d sources.\n", points.size());
    
    // samplers 
    std::random_device seed_generator;
    Sampler *main_sampler= nullptr;
    DebugSampler *debug_sampler= nullptr;
    
    size_t main_samples= size_t(options.width) * size_t(options.height) * size_t(options.samples_per_pixel);
    size_t pixel_samples= options.samples_per_pixel;
    int width= options.width;
    int height= options.height;
    
    if(options.debug_pixel_filename && !options.debug_export_samples)
    {
        width= options.debug_pixel_size * 2 +1;
        height= options.debug_pixel_size * 2 +1;
        main_samples= options.samples_per_pixel * width*height;
    }
    
    if(main_samples == 0)
    {
        main_samples= options.samples;
        pixel_samples= options.samples / (width * height);
    }
    
    assert(main_samples != 0);
    
    printf("using %lld samples: ", main_samples);
    printf("%dx%d %dspp\n", width, height, int(pixel_samples));
    
    //~ unsigned int seed= seed_generator();
    unsigned int seed= 0;
    
    if(std::string(options.sampler_name) == "samples")
        //~ main_sampler= new SamplerNd(options.samples_filename, options.dimensions, options.samples_remap_string, options.samples_pad_string, options.samples, /* seed */ seed);
        main_sampler= new SamplerNd(options.samples_filename, options.dimensions, nullptr, nullptr, options.samples, /* seed */ seed);  // remap / pad dimensiosn dans le sampler global
    
    // sampler sobol, sample_pack
    if(std::string(options.sampler_name) == "sobol")
        main_sampler= new SobolSampler(options.dimensions, main_samples, /* seed */ seed);
    
    // sampler sobol, sample_pack + owen scrambling
    if(std::string(options.sampler_name) == "owen")
        //~ main_sampler= new OwenSampler(options.dimensions, main_samples, /* seed */ seed);
        main_sampler= new OwenSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);
    
    // sampler sobol++ owen++ (par pixel)
    if(std::string(options.sampler_name) == "sobolpp")
        main_sampler= new SobolPPSampler(options.dimensions_remap, width, height, pixel_samples, options.seeds_filename, options.sobol_filename, options.owen_bits, /* seed */ seed);

    // sampler sobol++ owen (global)
    if(std::string(options.sampler_name) == "owenpp")
        main_sampler= new OwenPPSampler(options.dimensions_remap, main_samples, options.sobol_filename, /* seed */ seed);

    // sampler ZSampler hash
    if(std::string(options.sampler_name) == "zhash")
        main_sampler= new ZHashSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);
    // sampler ZSampler 
    if(std::string(options.sampler_name) == "z")
        main_sampler= new ZSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);

    // sampler Morton, re-impl de ZSampler 
    if(std::string(options.sampler_name) == "morton")
        main_sampler= new MortonSampler(options.dimensions, width, height, pixel_samples, options.sobol_filename, options.sobol_shift, /* seed */ seed);
    
    //~ // sampler ZPPSampler 
    //~ if(std::string(options.sampler_name) == "zpp")
        //~ main_sampler= new ZPPSampler(options.dimensions, width, height, pixel_samples, options.seeds_filename, options.sobol_filename, options.owen_bits, /* seed */ seed);
        
    // reference, siggraph 2019
    if(std::string(options.sampler_name) == "sampler19")
        main_sampler= new Sampler19(options.dimensions_remap, width, height, pixel_samples, /* seed */ seed);
    
    // reference, random + cr rotation
    if(std::string(options.sampler_name) == "random")
        main_sampler= new UniformSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);
    
    // reference, random
    if(std::string(options.sampler_name) == "rng")
        main_sampler= new RNGSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);
    
    // latin hyper square
    if(std::string(options.sampler_name) == "lhs")
        main_sampler= new LHSSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);
        
    // debug
    if(std::string(options.sampler_name) == "debug")
    {
        debug_sampler= new DebugSampler(options.dimensions, options.width, options.height, pixel_samples, seed);
        main_sampler= debug_sampler;
        
        options.debug_pixel_x= options.debug_pixel_x * options.width / 1024;
        options.debug_pixel_y= options.debug_pixel_y * options.height / 1024;
    }
    
    assert(main_sampler != nullptr);
    
    if(options.samples_remap_string)
        main_sampler->remap_dimensions(options.samples_remap_dimensions, options.samples_remap);
    else
        main_sampler->remap_dimensions(options.dimensions, options.dimensions_remap);
        
    //~ if(options.debug_pixel_filename || options.force_rotation)
    if(options.force_rotation)
        main_sampler->create_random_rotation(width, height, pixel_samples, options.rotation_scale);
    else if(options.force_null_rotation)
        main_sampler->force_null_rotation(width, height, pixel_samples);
    
    if(options.force_set_rotation)
        main_sampler->create_set_rotation(width, height, pixel_samples);
    
    if(options.force_tile)
        main_sampler->create_tile(options.tile_width, options.tile_height, width, height, pixel_samples);
    
#if 1
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
    
    {
        // echantillons nulls sur le domaine d'integration
        path_depth.resize(options.path_length+1, 0);
        last_depth.resize(options.width*options.height, 0);
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
    
    // initialise un sampler par run
    std::vector<Color> debug_pixel_runs;
    Framebuffer mean(image.width(), image.height());
    Framebuffer var(image.width(), image.height());
    
    for(int r= 0; r < options.runs; r++)
    {
        printf("%d threads.\n", nthreads);
        
        unsigned int s= seed;
        if(r > 0)
        {
            s= seed_generator();
            printf("run %d/%d seed %x\n", r, options.runs, s);
            // regenere les rotations cranley-patterson 
            main_sampler->seed(s);
        }
        
        // initialise un sampler par thread
        std::vector<Sampler *> samplers;
        for(int i= 0; i < nthreads; i++)
            samplers.push_back(main_sampler->clone(s));
        
        // naif, mais robuste...
        // \todo au lieu de stocker une image complete par thread, allouer un std::vector de paires (pixel, couleur)
        std::vector<Framebuffer *> framebuffers;
        for(int i= 0; i < nthreads; i++)
            framebuffers.push_back( new Framebuffer(options.width, options.height) );
        
        //~ const int next= 1;      // path tracing
        const int next= pixel_samples;      // rank / distance 
        size_t next_output= 0;
    #ifdef THREADS
        // decoupe l'espace d'iteration sur les samples
        {
            size_t samples= 0;
        #pragma omp parallel for schedule (dynamic, 4096)
            for(int64_t i= 0; i < int64_t(main_samples); i+= next)
            {
            #pragma omp atomic
                samples+= next;
                
                if(100*samples / main_samples != 100*(samples+next) / main_samples)
                {
                    printf("%lld/%d %lld%%...\r", samples / (pixel_samples * image.width()), image.height(), 100*samples / main_samples);
                    fflush(stdout);
                }
                
                // recupere le sampler associe au thread
                unsigned int tid= omp_get_thread_num();
                assert(tid < samplers.size());
                assert(tid < framebuffers.size());
                Sampler& sampler= *samplers[tid];
                Framebuffer& image= *framebuffers[tid];
                
                // indice du sample
                sampler.index(i);
                
    #else
        {
            // utilise uniquement le sampler principal
            Sampler& sampler= *samplers[0];
            Framebuffer& image= *framebuffers[0];
            
            // debug, pas de parallelisation
            for(size_t i= 0; i < sampler.size(); i++)
            {
                // indice du sample
                sampler.index(i);
    #endif
                
                // associe un pixel a l'indice
                int pi= i / pixel_samples;
                int px= pi % width;
                int py= pi / width;                
                Float x= Float(px) + Float(.5);
                Float y= Float(py) + Float(.5);
                
                if(options.sample_subpixel)
                {
                    sampler.dimension(0);       // force les dimensions utilisees pour generer le rayon primaire
                    x= sampler.sample1();
                    y= sampler.sample2();
                    
                    x= x * Float(width);
                    y= y * Float(height);
                    
                    if(options.debug_pixel_filename && !options.debug_export_samples)
                    {
                        x= x + options.debug_pixel_x - options.debug_pixel_size;
                        y= y + options.debug_pixel_y - options.debug_pixel_size;
                    }
                }
                
                if(options.debug_pixel_filename && options.debug_export_samples)
                {
                    px= std::floor(x);
                    py= std::floor(y);
                    if(px != options.debug_pixel_x || py != options.debug_pixel_y)
                        continue;
                        
                    // tres moche... enumere brutalement les samples associes au pixel 
                }
                
            #if 1
                if(options.force_pixel_center)
                {
                    // sampler global, sample selectionne le pixel, mais force le sample sur le centre du pixel / strate
                    px= std::floor(x);
                    py= std::floor(y);
                    x= px + Float(0.5);
                    y= py + Float(0.5);
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
                            {
                                // path tracing
                                //~ sample= sample + path(scene, time, /*depth */ options.path_length, points, sampler, ray) / lens_pdf;
                                sample= sample + path1_rank(scene, time, points, sampler, ray) / lens_pdf;
                            }
                        }
                        else
                        {
                            // pas de sources, ambient occlusion
                            sample= sample + ambient(scene, time, sampler, ambient_direction, ray) / lens_pdf;
                        }
                    }
                }
                
                // accumule par strate
                image.splat(x, y, Color(sample, 1));
            #endif
            
                //~ // enregistre les images intermediaires
                //~ if(options.force_output && !options.debug_pixel_filename && i == next_output)
                //~ {
                    //~ std::chrono::high_resolution_clock::time_point cpu_stop= std::chrono::high_resolution_clock::now();
                    //~ int cpu_time= std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
                    //~ printf("%04d/%d samples, cpu  %ds %03dms\n", int(s)+1, options.samples_per_pixel, (int) (cpu_time / 1000), (int) (cpu_time % 1000));

                    //~ Framebuffer tmp(image.width(), image.height());
                    //~ for(int tid= 0; tid < int(framebuffers.size()); tid++)
                        //~ tmp.add(framebuffers[tid]);
                    
                    //~ write_image_hdr(tmp.flatten(), Format("%s-%d-%04d.hdr", options.output_prefix, options.samples_per_pixel, int(i) +1));
                    //~ write_image(tmp.flatten(2.2), Format("%s-%d-%04d.png", options.output_prefix, options.samples_per_pixel, int(i) +1));
                    
                    //~ next_output= (i+1)*2 -1;
                //~ }
                
                //~ if(options.force_output && options.debug_pixel_filename && i == next_output)
                //~ {
                    //~ Color pixel= image(options.debug_pixel_x, options.debug_pixel_y);
                    //~ printf("  %lu %.15lf %.15lf %.15lf\n", i+1, pixel.r, pixel.g, pixel.b);
                    
                    //~ next_output= (i+1)*2 -1;
                //~ }
            }   // for(sample...)
        }
        printf("\ndone.\n");
        
        //
        for(int tid= 0; tid < int(framebuffers.size()); tid++)
            image.add(framebuffers[tid]);
        
        if(options.runs > 1)
        {
            if(!options.debug_pixel_filename)   // pas la peine de sauver l'image pour quelques pixels...
                write_image_pfm(image.flatten(), Format("%s-%05lu-run-%03d.pfm", options.output_prefix, pixel_samples, r));
            
            // accumule les runs
            if(options.debug_reference)
            {
                // erreur moyenne
                for(unsigned i= 0; i < image.size(); i++)
                {
                    if(image.pixels[i].grey() == 0)
                        continue;
                    
                    Color32 ref= debug_reference(i);
                    Color error= Color(ref.r, ref.g, ref.b) - image.pixels[i];
                    error= Color(std::abs(error.r), std::abs(error.g), std::abs(error.b));
                    mean.pixels[i]= mean.pixels[i] + error / Float(options.runs);
                }
            }
            else
            {
                // moyenne des runs
                for(unsigned i= 0; i < image.size(); i++)
                {
                    mean.pixels[i]= mean.pixels[i] + image.pixels[i] / Float(options.runs);
                    var.pixels[i]= var.pixels[i] + image.pixels[i]*image.pixels[i] / Float(options.runs);
                }
            }
            
            if(options.debug_pixel_filename)
            {
                Color pixel= image(options.debug_pixel_x, options.debug_pixel_y);
                debug_pixel_runs.push_back( pixel );
            }
            
            // remise a zero de l'image pour le prochain run
            if(r < options.runs -1)
                image.clear();  // sauf le dernier...
        }
        
        // prepare les samplers pour le prochain run
        for(int tid= 0; tid < int(samplers.size()); tid++)
            delete samplers[tid];
        
        for(int tid= 0; tid < int(framebuffers.size()); tid++)
            delete framebuffers[tid];
    }   // for(run...)
    
    // affiche quelques stats
    std::chrono::high_resolution_clock::time_point cpu_stop= std::chrono::high_resolution_clock::now();
    int cpu_time= std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
    printf("cpu  %ds %03dms\n", (int) (cpu_time / 1000), (int) (cpu_time % 1000));
    printf("%fMsamples/s\n", (image.width() * image.height() * options.samples_per_pixel / Float(1000000.0)) / (Float(cpu_time) / Float(1000.0)));
    
    // enregistre les images
    //~ write_image_hdr(image.flatten(), Format("%s-%05lu.hdr", options.output_prefix, pixel_samples));
    write_image_pfm(image.flatten(), Format("%s-%05lu.pfm", options.output_prefix, pixel_samples));
    //~ write_image_exr(image.flatten(), Format("%s-%05lu.exr", options.output_prefix, pixel_samples));

    //~ printf("writing map image '%s'...\n", (const char *) Format("%s-%05lu-map.hdr", options.output_prefix, pixel_samples));
    //~ write_image_hdr(image.flatten_weights(), Format("%s-%05lu-map.hdr", options.output_prefix, pixel_samples));
    
    //~ printf("writing zeroes image '%s'...\n", (const char *) Format("%s-%05lu-zeroes.hdr", options.output_prefix, pixel_samples));
    //~ write_image_hdr(image.flatten_zeroes(), Format("%s-%05lu-zeroes.hdr", options.output_prefix, pixel_samples));
    
    //~ printf("writing output image '%s'...\n", (const char *) Format("%s-%05lu.png", options.output_prefix, main_samples));
    //~ write_image(image.flatten(Float(2.2)), Format("%s-%05lu.png", options.output_prefix, main_samples));
    
//    printf("writing output image '%s'...\n", options.output_filename);
//    write_image(image.flatten(2.2), options.output_filename);
    
    {
        // chemins nulls / non nulls
        size_t total= 0;
        for(int d= 0; d < int(path_depth.size()); d++)
        {
            total+= path_depth[d];
            float paths= float(path_depth[d]) / float(options.width * options.height);
            printf("depth %d: samples %lu, samples per pixel %f %.2f%%\n", d, path_depth[d], paths, 100 * paths / pixel_samples);
        }
    }
    
    if(options.debug_pixel_filename && options.runs > 1)
    {
        Color mean;
        for(auto& color : debug_pixel_runs)
            mean= mean + color / Float(debug_pixel_runs.size());
        
        Color var;
        for(auto& color : debug_pixel_runs)
            var= var + ((color - mean) * (color - mean)) / Float(debug_pixel_runs.size());
        
        FILE *out= fopen(options.debug_pixel_filename, "wt");
        
        for(auto& color : debug_pixel_runs)
            fprintf(out, "%llu %.15lf %.15lf %.15lf\n", pixel_samples, double(color.r), double(color.g), double(color.b));
        
        fprintf(out, "# mean %llu %.15lf %.15lf %.15lf\n", debug_pixel_runs.size(), double(mean.r), double(mean.g), double(mean.b));
        fprintf(out, "# var %llu %.15lf %.15lf %.15lf\n", debug_pixel_runs.size(), double(var.r), double(var.g), double(var.b));
        fclose(out);
        
        printf("writing debug pixel '%s'... %.15lf %.15lf %.15lf\n", options.debug_pixel_filename, mean.r, mean.g, mean.b);
    }
    if(options.runs > 1)
    {
        // accumule les runs
        for(unsigned i= 0; i < image.size(); i++)
            var.pixels[i]= var.pixels[i] - mean.pixels[i]*mean.pixels[i];
            
        write_image_pfm(mean.flatten(), Format("%smean-%05lu.pfm", options.output_prefix, pixel_samples));
        //~ write_image_pfm(var.flatten(), Format("%svar-%05lu.pfm", options.output_prefix, pixel_samples));
    }
    
    main_sampler->release();
    delete main_sampler;
#endif

    return 0;
}
