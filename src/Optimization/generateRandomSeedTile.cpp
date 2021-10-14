#include <CLI11.hpp>
#include <SeedTileHelper.hpp>
#include <random>
#include <assert.h>

// Generate a random seed tile and export it
// (+ debug to check the save/load a map)

int main(int argc, char **argv)
{
  CLI::App app { "Generate a random seed tile."};
  unsigned int size = 128;
  app.add_option("-s,--size", size, "Size of the image grid (default 128)");
  uint32_t seed = 12345;
  app.add_option("--seed", seed, "Seed (def 12345)");
  uint32_t dimension = 2;
  app.add_option("-d,--dimension", dimension, "Dimension of the seedmap (def 2)");
  std::string output_fname;
  app.add_option("-o,--output", output_fname, "Output filename .")->required();
  CLI11_PARSE(app, argc, argv)
    
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);  // uniform distribution of integers between 0 and 2^32-1
  
  SeedTile seedmap(size,dimension);
  
  for(auto i=0 ; i < size; ++i)
  for(auto j=0; j < size; ++j)
  for(auto d = 0; d < dimension; ++d)
  {
    seedmap(i,j, d) = unif32bits(gen);
  }
  std::cout<<"Seed map generated."<<std::endl;
  seedmap.saveTile(output_fname);
  
  //For debug, just loading and checking the values
  SeedTile seedmap2(size, dimension);
  seedmap2.loadTile(output_fname);
  for(auto i=0; i < size; ++i)
    for(auto j=0; j < size; ++j)
    for(auto d = 0; d < dimension; ++d)
     assert( seedmap(i,j,d) == seedmap2(i,j,d));
  return 0;
}
