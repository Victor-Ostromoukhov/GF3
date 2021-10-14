
#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>



struct Options
{
    //~ const char *file_filename= "scenes/cornell3_bars.obj";     // objet simple sans description de scene
    //~ const char *scene_filename= nullptr;    // fichier scene
    //~ const char *orbiter_filename= "scenes/cornell3_orbiter.txt";
    //~ const char *lens_filename= "scenes/cornell_lens.txt";
    
    const char *file_filename= nullptr;     // objet simple sans description de scene
    const char *scene_filename= nullptr;    // fichier scene
    const char *orbiter_filename= nullptr;
    const char *lens_filename= nullptr;
    
    const char *samples_filename= "<no file>";
    const char *samples_remap_string= nullptr;
    std::vector<int> samples_remap= std::vector<int>();
    int samples_remap_dimensions= 0;
    
    const char *dimensions_remap_string= nullptr;
    std::vector<int> dimensions_remap= std::vector<int>();
    int dimensions= 7;
    
    const char *output_filename= "render.png";
    const char *output_prefix= "render";
    
    const char *sampler_name= "sobol";
    const char *sobol_filename= nullptr;
    const char *seeds_filename= nullptr;
    int owen_bits= 32;
    int sobol_shift= 0;
    
    int width= 1024;
    int height= 1024;
    int samples_per_pixel= 1;
    size_t samples= 0;
    int path_length= 0;
    
    int tile_width= 0;
    int tile_height= 0;
    
    const char *debug_pixel_filename= nullptr;
    const char *debug_reference= nullptr;
    int debug_pixel_runs= 1;
    int debug_pixel_x= -1;
    int debug_pixel_y= -1;
    int debug_pixel_size= 3;
    bool debug_export_samples= false;
    
    int runs= 1;
    int threads= 0;
    
    int rotation_scale= 0;
    
    bool force_pixel_center= false;
    bool force_rotation= false;
    bool force_null_rotation= false;
    bool force_set_rotation= false;
    bool force_tile= false;
    
    bool force_output= false;
    bool force_last_bounce= false;
    
    bool force_diffuse_material= true;
    bool force_direct_lighting= true;
    
    bool remap_bounce= false;
    bool remap_bounce3d= false;
    bool remap_bounce2d= false;
    bool source_2d= false;
    bool brdf_2d= false;
    bool remap1= false;
    bool source_quads= false;
    bool source_quad1= false;
    
    bool no_shift= false;
    bool sample_subpixel= false;
    bool sample_lens= false;
    bool lens_force_square= true;
    bool sample_time= false;
    
    bool force_ambient_occlusion= false;
    bool ambient_use_disk_cosine= false;
    bool ambient_use_disk_uniform= true;
    bool ambient_use_uniform= false;
    bool ambient_use_fibonacci= false;
    
    Options() = default;
    Options( const Options& ) = delete;
    Options& operator= ( const Options& ) = default;
    
    Options( const char **begin, const char **end ) : options(begin, end)
    {
        if(options.size() > 1)
        {
            //~ file_filename= value("--file", "scenes/cornell3_bars.obj");     // pas de valeur par defaut
            file_filename= value("--file", nullptr);
            scene_filename= value("--scene", nullptr);
            
            //~ orbiter_filename= value("--camera", "scenes/cornell3_orbiter.txt");     // pas de valeur par defaut
            orbiter_filename= value("--camera", nullptr);
            if(find_value("--lens"))
            {
                //~ lens_filename= value("--lens", "scenes/cornell_lens.txt");  // pas de valeur par defaut
                lens_filename= value("--lens", nullptr);
                sample_lens= true;
            }
            sample_lens= flag("--sample_lens");
            
            output_filename= value("-o", "render.png");
            output_prefix= value("--prefix", "render6d");
            
            width= value_int("-w", "1024");
            height= value_int("-h", "1024");
            samples_per_pixel= value_int("-s", "0");
            sample_subpixel= flag("--sample_subpixel");
            path_length= value_int("--path_length", "1");
            sample_time= flag("--sample_time");
            
            if(find_value("--samples"))
            {
                samples_filename= value("--samples", "samples.txt");
                sampler_name= "samples";
                dimensions= value_int("-d", "7");
                samples= value_long("-n", "0");
                // samples_per_pixel= 0;
                sample_subpixel= true;      // sampler global, distribution arbitraire sur le plan image
            }
            else
            {
                sampler_name= value("--sampler", "random");
                dimensions= value_int("-d", "7");
                samples= value_long("-n", "0");
            }
            
            samples_remap_string= value("--remap_dims", nullptr);
            samples_remap= parse_dimensions(samples_remap_string);
            if(samples_remap.size() > 0)
                samples_remap_dimensions= samples_remap.size();
            
            dimensions_remap_string= value("--dims", nullptr);
            dimensions_remap= parse_dimensions(dimensions_remap_string);
            if(dimensions_remap.size() > 0)
                dimensions= dimensions_remap.size();
            
            sobol_filename= value("--sobol_table", nullptr);
            owen_bits= value_int("--owen_bits", "32");
            sobol_shift= value_int("--sobol_shift", "0");
            seeds_filename= value("--seeds", nullptr);
            //~ no_shift= flag("--no_shift");
            
            // force --sample_subpixel pour sobol, sampler global, distribution arbitraire sur le plan image
            if(std::string(sampler_name) == "sobol") 
                sample_subpixel= true;
            
            force_pixel_center= value_bool("--force_pixel_center", force_pixel_center ? "true" : "false");
            
            force_ambient_occlusion= flag("--force_ao");
            ambient_use_disk_cosine= flag("--ao_disk_cosine");
            ambient_use_disk_uniform= flag("--ao_disk_uniform");
            ambient_use_uniform= flag("--ao_uniform");
            
            force_diffuse_material= flag("--force_diffuse");
            force_direct_lighting= flag("--force_direct");
            force_last_bounce= flag("--force_last_bounce");
            
        #ifdef REMAP_BOUNCE 
            remap_bounce= true;
        #endif
        #ifdef BRDF_2D 
            brdf_2d= true;
        #endif
        #ifdef SOURCE_2D
            source_2d= true;
        #endif
        #ifdef SOURCE_QUADS
            source_quads= true;
        #endif
        #ifdef SOURCE_QUAD1
            source_quad1= true;
        #endif            
            remap_bounce= value_bool("--remap_bounce", remap_bounce ? "true" : "false");
            source_2d= value_bool("--source_2d", source_2d ? "true" : "false");
            brdf_2d= value_bool("--brdf_2d", brdf_2d ? "true" : "false");
            source_quads= value_bool("--source_quads", source_quads ? "true" : "false");
            source_quad1= value_bool("--source_quad1", source_quad1 ? "true" : "false");
            
            // sanity check remapping options
            if(remap_bounce && source_2d && brdf_2d)
                remap_bounce2d= true;
            else if(remap_bounce)
            {
                remap_bounce3d= true;
                brdf_2d= false,
                source_2d= false;
            }
            
            force_output= flag("--force_output");

            runs= value_int("--runs", "1");
            threads= value_int("--threads", "0");
            
            if(find_value("--debug_pixel", 4))
            {
                int index= find("--debug_pixel");
                assert(index != -1);
                
                debug_pixel_filename= options[index +1];
                if(sscanf(options[index+2], "%d", &debug_pixel_x) != 1)
                    printf("invalid parameter '%s'\n", options[index+2]);
                if(sscanf(options[index+3], "%d", &debug_pixel_y) != 1)
                    printf("invalid parameter '%s'\n", options[index+3]);
                if(sscanf(options[index+4], "%d", &debug_pixel_runs) != 1)
                    printf("invalid parameter '%s'\n", options[index+4]);
                
                options.erase(options.begin() + index, options.begin() + index +5);
                
                //
                runs= debug_pixel_runs;
                
                // 
                debug_export_samples= flag("--export_samples");
                debug_pixel_size= value_int("--export_size", "0");
                debug_reference= value("--debug_reference", nullptr);
            }
            
            force_rotation= flag("--force_rotation");
            force_null_rotation= flag("--force_null_rotation");
            force_set_rotation= flag("--force_set_rotation");
            
            rotation_scale= value_int("--rotation_scale", "0");
            
            if(find_value("--force_tile", 2))
            {
                int index= find("--force_tile");
                assert(index != -1);
                
                force_tile= true;
                if(sscanf(options[index +1], "%d", &tile_width) != 1)
                    printf("invalid parameter '%s'\n", options[index +1]);
                if(sscanf(options[index +2], "%d", &tile_height) != 1)
                    printf("invalid parameter '%s'\n", options[index +2]);
                
                options.erase(options.begin() + index, options.begin() + index +3);
            }
        }
        
        if(options.size() == 1 || flag("--help"))
        {
            printf("configuration:\n");
            printf("  --threads %d, openmp max threads\n", threads);
            
            if(scene_filename)
                printf("  --scene '%s'\n", scene_filename);
            else
                printf("  --file '%s'\n", file_filename);
            printf("  --camera '%s'\n", orbiter_filename);
            printf("  --lens '%s'\n", lens_filename);
            
            printf("  --prefix '%s'\n", output_prefix);
            printf("  -o '%s'\n", output_filename);
            printf("  -w %d, width\n", width); 
            printf("  -h %d, height\n", height); 
            printf("  -s %d, samples per pixel\n", samples_per_pixel); 
            printf("  -n %llu, samples (total per image)\n", samples);
            printf("  --sampler '%s'\n", sampler_name);
            printf("  --samples '%s'\n", samples_filename);
            printf("  --sobol_table '%s'\n", sobol_filename ? sobol_filename : "<default file>");
            printf("  --sobol_shift '%d', 0 to use no shift\n", sobol_shift);
            printf("  --owen_bits %d, 0 to use log2(nspp) bits\n", owen_bits);
            printf("  -d %d, dimensions\n", dimensions);
            printf("  --dims \"dimensions\" '%s'\n", dimensions_remap_string);
            if(dimensions_remap_string)
            {
                printf("    explicit %d dimensions: ", dimensions);
                for(int i= 0; i < dimensions; i++) 
                    printf("(index %d= dim %d) ", i, dimensions_remap[i]);
                printf("\n");
            }
            else
            {
                printf("    ");
                for(int i= 0; i < dimensions; i++) 
                    printf("(index %d= dim %d) ", i, dimensions_remap[i]);
                printf("\n");
            }
            
            printf("  --remap_dims \"dimensions\" '%s'\n", samples_remap_string);
            if(samples_remap_string)
            {
                printf("    remapping %d dimensions: ", samples_remap_dimensions);
                for(int i= 0; i < samples_remap_dimensions; i++) 
                    printf("(sample[%d]= index %d= dim %d) ", i, samples_remap[i], (samples_remap[i] != -1) ? dimensions_remap[samples_remap[i]] : -1);
                printf("\n");
            }
            
            printf("  --sample_subpixel '%s'\n", sample_subpixel ? "true" : "false");
            printf("  --force_pixel_center '%s'\n", force_pixel_center ? "true" : "false");
            printf("  --sample_time '%s'\n", sample_time ? "true" : "false");
            printf("  --sample_lens '%s', lens file '%s'\n", sample_lens ? "true" : "false", lens_filename);
            printf("  --path_length %d\n", path_length);
            printf("  --force_last_bounce '%s'\n", force_last_bounce ? "true" : "false");
            printf("  --force_direct '%s', render direct lighting\n", force_direct_lighting ? "true" : "false");
            
            printf("  --remap_bounce '%s', use same random numbers for sampling light sources and sampling brdf\n", remap_bounce ? "true" : "false");
            printf("  --source_2d '%s', use 2 random numbers + renormalization to sample light sources\n", source_2d ? "true" : "false");
            printf("  --brdf_2d '%s', use 2 random numbers + renormalization to sample brdfs\n", brdf_2d ? "true" : "false");
            printf("  --source_quads '%s', use quad light sources\n", source_quads ? "true" : "false");
            printf("  --source_quad1 '%s', use a single quad light source\n", source_quad1 ? "true" : "false");
            
            printf("  --force_diffuse '%s', use diffuse materials\n", force_diffuse_material ? "true" : "false");
            
            printf("  --force_ao '%s', render ambient occlusion\n", force_ambient_occlusion ? "true" : "false");
            printf("  --ao_disk_cosine '%s', use concentric disk cosine directions\n", ambient_use_disk_cosine ? "true" : "false");
            printf("  --ao_disk_uniform '%s', use concentric disk uniform directions\n", ambient_use_disk_uniform ? "true" : "false");
            printf("  --ao_uniform '%s', use uniform directions\n", ambient_use_uniform ? "true" : "false");
            
            printf("  --force_output '%s', write temp results\n", force_output ? "true" : "false");
            
            printf("  --debug_pixel '%s', output '%s', x %d y %d, runs %d\n", debug_pixel_filename ? "true" : "false", debug_pixel_filename ? debug_pixel_filename  : "<no file>", 
                debug_pixel_x, debug_pixel_y, debug_pixel_runs);
            
            printf("  --runs %d, render several images with different seeds\n", runs);
            
            printf("  --force_rotation '%s', use same sequence per pixel + cranley-patterson rotation\n", force_rotation ? "true" : "false");
            printf("  --rotation_scale %d, use micro scale cranley-patterson rotation\n", rotation_scale);
            printf("  --force_null_rotation '%s', use same sequence per pixel\n", force_null_rotation ? "true" : "false");
            printf("  --force_set_rotation '%s', use one sequence per pixel, several sets\n", force_set_rotation ? "true" : "false");
            printf("  --force_tile %d %d, use one sequence per tile\n", tile_width, tile_height);
            
            printf("\n");
            
            printf("command line:\n");
            for(const char **a= begin; a != end; a++)
                printf("%s ", *a);
            printf("\n");
            
            printf("unprocessed options:\n");
            for(int i= 0; i < int(options.size()); i++)
                printf("%s ", options[i]);
            printf("\n\n");
        }
    }

    
    std::vector<int> parse_dimensions( const char *remap_string )
    {
        std::vector<int> remap;
        int ndim= dimensions;
        
        // dimension padding
        if(remap_string)
        {
            int offset= 0;
            int next= 0;
            ndim= 0;
            if(sscanf(remap_string + offset, "%d %n", &ndim, &next) == 1)
            {
                for(int i= 0; i < ndim; i++)
                    remap.push_back(-1);
                    
                for(int d= 0; d < ndim; d++)
                {
                    int rd= 0;
                    int rdim= 0;
                    offset= offset + next;
                    if(sscanf(remap_string + offset, "%d:%d %n", &rdim, &rd, &next) != 2)
                        break;
                    
                    //~ if(rd < 0 || rd >= ndim)
                    //~ {
                        //~ printf("[error] can't map data[%d] to sample[%d], %dd data\n", rd, rdim, ndim);
                        //~ remap.clear();
                        //~ break;
                    //~ }
                    if(rdim < 0 || rdim >= ndim)
                    {
                        printf("[error] can't map index %d to dim %d, %dd index\n", rd, rdim, ndim);
                        remap.clear();
                        break;
                    }
                    
                    remap[rdim]= rd;
                }
            }
        }
        else
        {
            assert(dimensions > 0);
            // identity mapping
            for(int i= 0; i < dimensions; i++)
                remap.push_back(i);
        }
        
        if(remap_string && remap.empty())
        {
            printf("[error] remapping dimensions: '%s'...\n", remap_string);
            exit(1);
        }
        
        printf("option '%s'\n", remap_string);
        printf("  %d dimensions: ", ndim);
        for(int i= 0; i < ndim; i++) 
            printf("(index %d= dim %d) ", i, remap[i]);
        printf("\n");
        
        return remap;
    }
    
    // recherche des options
    int find( const char *name )
    {
        for(int i= 0; i < int(options.size()); i++)
            if(std::string(options[i]) == std::string(name))
                return i;
        return -1;
    }

    // cherche un flag, renvoie un bool
    bool flag( const char *name )
    {
        int index= find(name);
        if(index != -1)
            options.erase(options.begin() + index, options.begin() + index +1);
        return (index != -1);
    }

    bool find_value( const char *name, const int n= 1 )
    {
        int index= find(name);
        return (index >= 0 && index + n < int(options.size()));
    }
    
    // cherche une option et renvoie sa valeur ou la valeur par defaut
    const char *value( const char *name, const char *default_value )
    {
        int index= find(name);
        if(index < 0 || index + 1 >= int(options.size()))
            return default_value;
        
        const char *str= options[index +1];
        options.erase(options.begin() + index, options.begin() + index +2);
        return str;
    }

    // cherche une option, renvoie un int
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
    
    // cherche une option, renvoie un unsigned long int 
    size_t value_long( const char *name, const char *default_value )
    {
        const char *str= value(name, default_value);
        
        unsigned long int v= 0;
        if(sscanf(str, "%lu", &v) == 1)
            return v;
        
        // erreur parametre invalide.
        printf("invalid int parameter: %s = '%s'\n", name, str);
        
        exit(1);
        return 0;
    }

    // cherche une option, renvoie un bool
    bool value_bool( const char *name, const char *default_value )
    {
        const char *str= value(name, default_value);
        if(std::string(str) == "true" || std::string(str) == "on" || std::string(str) == "1")
            return true;
        if(std::string(str) == "false" || std::string(str) == "off" || std::string(str) == "0")
            return false;
        
        // erreur parametre invalide.
        printf("invalid bool parameter: %s = '%s'\n", name, str);
        
        exit(1);
        return false;
    }
    
protected:
    std::vector<const char *> options = std::vector<const char *>();
};

extern Options options;

#endif
    
    
