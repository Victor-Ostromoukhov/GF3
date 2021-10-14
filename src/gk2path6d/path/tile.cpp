
#include "sampler.h"
#include "sampler_sobol.h"
#include "sampler_sobolpp.h"
#include "sampler_zsampler.h"
#include "sampler_morton.h"

#include "options.h"

Options options;

int main( int argc, char *argv[] )
{
    options= Options((const char **) argv, (const char **) argv + argc);

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
    printf("using %lld samples: ", main_samples);
    printf("%dx%d %dspp\n", width, height, int(pixel_samples));
    
    unsigned int seed= 0;
    Sampler *main_sampler= nullptr;
    
    // sampler sobol, sample_pack
    if(std::string(options.sampler_name) == "sobol")
        main_sampler= new SobolSampler(options.dimensions, main_samples, /* seed */ seed);
    
    if(std::string(options.sampler_name) == "owen")
        //~ main_sampler= new OwenSampler(options.dimensions, main_samples, /* seed */ seed);
        main_sampler= new OwenSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);
    
    // sampler sobol++ owen++ (par pixel)
    if(std::string(options.sampler_name) == "sobolpp")
        main_sampler= new SobolPPSampler(options.dimensions_remap, width, height, pixel_samples, options.seeds_filename, options.sobol_filename, options.owen_bits, /* seed */ seed);

    // sampler sobol++ owen (global)
    if(std::string(options.sampler_name) == "owenpp")
        main_sampler= new OwenPPSampler(options.dimensions_remap, main_samples, options.sobol_filename, /* seed */ seed);

    // sampler ZSampler 
    if(std::string(options.sampler_name) == "z")
        main_sampler= new ZSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);

    // sampler Morton, re-impl de ZSampler 
    if(std::string(options.sampler_name) == "morton")
        main_sampler= new MortonSampler(options.dimensions, width, height, pixel_samples, options.sobol_filename, options.sobol_shift, /* seed */ seed);

    // sampler Morton01, re-impl de ZSampler 
    if(std::string(options.sampler_name) == "morton01")
        main_sampler= new MortonSampler01(options.dimensions, width, height, pixel_samples, options.sobol_filename, options.sobol_shift, /* seed */ seed);
    
    // reference, random + cr rotation
    if(std::string(options.sampler_name) == "random")
        main_sampler= new UniformSampler(options.dimensions, width, height, pixel_samples, /* seed */ seed);

    assert(main_sampler);
    
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
    
    FILE *out= fopen(options.output_filename, "wt");
    if(!out)
    {
        printf("[error] creating output file '%s'...\n", options.output_filename);
        exit(1);
    }
    
    // genere les samples pour la tuile 
    int ndim= main_sampler->dimensions();
    for(size_t i= 0; i < main_samples; i++)
    {
        main_sampler->index(i);
        for(int d= 0; d < ndim; d++)
            fprintf(out, "%.10f ", main_sampler->sample1());
        fprintf(out, "\n");
    }
    
    fclose(out);
    printf("writing tile '%s'...\n", options.output_filename);
    return 0;
}
