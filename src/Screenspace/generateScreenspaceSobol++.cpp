#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>

#include <CLI11.hpp>
#include <SeedTileHelper.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

#include "../Integration/Integration.h"


using namespace std;

struct t_int_point2D {
	uint32_t x, y;
};

inline float get_pixel_sample(  const std::vector<SobolGenerator1D<uint32_t> >& sobols,
                              const SeedTile &seedMap,
                              const uint32_t i,
                              const uint32_t j,
                              const uint32_t spp,
                              const uint32_t dim)
{
  return getOwenPlus1D_with_seed(sobols, dim, spp, seedMap(i % seedMap.size,j % seedMap.size,dim)); 	// getOwenPlus1D_with_seed() is in src/Samplers/SobolGenerator1D.h
}

uint32_t seed = 13374269;
std::mt19937_64 gen(seed);
std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);  // uniform distribution of integers between 0 and 2^32-1


int main(int argc, char **argv)
{
  CLI::App app { "genScreenspace nD using sobol++. Generate samples per pixel in nD."};
  uint size = 8;
  uint spp = 4;
  uint dim = 2;
  app.add_option("-s,--size", size, "Size of the image grid (default 8)");
  app.add_option("-n,--spp", spp, "Number of spp (default 4)");
  app.add_option("-d,--dim",dim, "Dimension (defaukt 2)");
  std::string output_fname;
  app.add_option("-o,--output", output_fname, "Output filename for the samples (ascii .Dat).")->required();
  std::string dir_vectors_fname; // = "../../../data/sobol_init_tab.dat";
  app.add_option("--dirs", dir_vectors_fname, "File name of the Sobol intialization table (e.g. ../../../data/sobol_init_tab.dat)")->required()->check(CLI::ExistingFile);
  std::string seed_map_fname;
  app.add_option("-t,--tile", seed_map_fname, "Filename of the seed map.")->required()->check(CLI::ExistingFile);
  std::vector<int> indices;
  app.add_option("--indices", indices, "Sobol index per dimension NOT USED YET (e.g. --indices 5 6");
  bool exportImages=false;
  app.add_flag("-x", exportImages, "Export first samples (1spp) as images (one per dimension)");
  bool skipPixelCoord=false;
  app.add_flag("--skipPixelCoord", skipPixelCoord, "Skip pixel coord when exporting");
  bool performIntegration=false;
  app.add_flag("--performIntegration", performIntegration, "Perform per pixel integration and export of the HDR");
  CLI11_PARSE(app, argc, argv)
  
  std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
  loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file and fill appropriate structures
  
  //TileMap
  SeedTile seedMap;
  seedMap.loadTile(seed_map_fname);
  if (seedMap.dimension < dim)
  {
    std::cout<<"[ERROR] the Seed map has dimension "<<seedMap.dimension<<" but you ask samples of dimension "<<dim<<std::endl;
    exit(1);
  }
  if (seedMap.dimension != dim)
    std::cerr<<"[WARNING] the Seed map has dimension "<<seedMap.dimension<<" but you ask samples of dimension "<<dim<<". I will only used the first dimensions of the seedMap"<<std::endl;
  if (seedMap.size < size)
    std::cerr<<"[WARNING] the tile map size is inferior to the requested size, I will do modulo "<< seedMap.size<<std::endl;
  
  std::ofstream ofs(output_fname, std::ios_base::out);
  
  //Generate the points
  for(auto i=0; i < size; ++i)
  {
    for(auto j=0; j < size; ++j)
    {
      for(auto sampleid=0; sampleid< spp ; sampleid++)
      {
        if (!skipPixelCoord)
          ofs<<i<<" "<<j<<" ";
        for(auto d= 0; d < dim; ++d)
           ofs << get_pixel_sample(sobols,seedMap, i,j, sampleid, d) <<" ";
        ofs<<std::endl;
      }
    }
  }
  ofs.close();
  
  if (exportImages)
  {
    std::cout<<"Exporting coordinates"<<std::endl;
    std::vector<float> result(size*size);
    for(auto d=0; d < dim; ++d)
    {
      std::string filename="output-"+std::to_string(d)+".hdr";
      for(auto i=0; i < size; ++i)
      {
        for(auto j=0; j < size; ++j)
        {
          result[i+j*size] = get_pixel_sample(sobols,seedMap, i,j, 0, d);
        }
      }
      stbi_write_hdr(filename.c_str(), size,size,1, result.data());
    }
  }
  
  if (performIntegration)
  {
    //Generate the points
    std::vector<float> result(size*size);
    for(auto i=0; i < size; ++i)
    {
      for(auto j=0; j < size; ++j)
      {
        std::vector< VecXDynamic > points;
        for(auto s = 0; s < spp ; ++s)
        {
          VecXDynamic sample(dim);
          for(auto d= 0; d < dim; ++d)
          sample[ d ] =  get_pixel_sample(sobols,seedMap, i,j, s, d) ;
          points.push_back( sample );
        }

          double v =  calculate_mse(points,2,2048);
          result[i+j*size] =  v ; //Gaussr
        }
      }
    stbi_write_hdr("export-MSE.hdr", size,size,1, result.data());
  }
  
  return 0;
}
