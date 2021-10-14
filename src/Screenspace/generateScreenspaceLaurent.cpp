#include <iostream>
#include <functional>

#include <CLI11.hpp>
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#ifdef _MSC_VER 
typedef unsigned int uint;
#endif

int main(int argc, char **argv)
{
  CLI::App app { "genScreenspace nD. Generate samples per pixel in nD." };
  uint size = 8;
  uint spp = 4;
  uint dim = 2;
  uint mask = 256;
  app.add_option("-s,--size", size, "Size of the image grid (default 8)");
  app.add_option("-n,--spp", spp, "Number of spp (default 4)");
  app.add_option("-d,--dim",dim, "Dimension (defaukt 2)");
  app.add_option("-m,--mask", mask, "Number of spp of the optimized mask to use (>=spp) (default=256)");
  bool exportImages=false;
  app.add_flag("-x", exportImages, "Export first samples as images");
  CLI11_PARSE(app, argc, argv)

  std::function<double(int,int,int,int)> sampler;
  switch (mask)
  {
    case   1: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp;   break;
    case   2: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp;   break;
    case   4: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp;   break;
    case   8: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp;   break;
    case  16: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp;  break;
    case  32: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp;  break;
    case  64: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp;  break;
    case 128: sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp; break;
    default:  sampler = samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp; break;
  }
  
  for(auto i=0; i < size; ++i)
  for(auto j=0; j < size; ++j)
  for(auto sampleid=0; sampleid< spp ; sampleid++)
  {
    std::cout<<i<<" "<<j<<" ";
    for(auto d= 0; d < dim; ++d)
      std::cout << sampler(i,j, sampleid, d) <<" ";
    std::cout<<std::endl;
  }
  
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
          result[i+j*size] = sampler(i,j, 0, d);
        }
      }
      stbi_write_hdr(filename.c_str(), size,size,1, result.data());
    }
  }
  
  return 0;
}
