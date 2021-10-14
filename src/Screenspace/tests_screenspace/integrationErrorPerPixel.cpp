#include <iostream>
#include <random>

#include <CLI11.hpp>  
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp.hpp"
#include "../samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp.hpp"

#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

#include "../Integration/Gaussians.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../includes/stb_image_write.h"

int main(int argc, char **argv)
{
  CLI::App app { "integrationError per pixel 2D" };
  uint size = 8;
  uint spp = 4;
  uint mask = 256;
  uint dimension = 2;
  std::string filename;
  app.add_option("-d,--dim", dimension, "Dimension of the integration problem (default 2)",true);
  app.add_option("-s,--size", size, "Size of the image grid (default 8)",true);
  app.add_option("-n,--spp", spp, "Number of spp (default 4)",true);
  app.add_option("-m,--mask", mask, "Number of spp of the optimized mask to use (>=spp) (default=256)",true);
  app.add_option("-o,--output", filename, "Output filename of the HDR containing the integration errors")->required();
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
  
  std::cout<<"=============================================="<<std::endl;
  std::cout<<"Computing pex pixel integration:"<<std::endl;
  std::cout<<app.config_to_str(true,true);
  std::cout<<"=============================================="<<std::endl;

  std::vector<VecX<2>> samples(spp);
  std::vector<float> result(size*size);
  double totalmse=0.0;
  
  for(auto i=0; i < size; ++i)
  for(auto j=0; j < size; ++j)
  {
    //samples (2D only at this point)
    for(auto sampleid=0; sampleid< spp ; sampleid++)
      samples[sampleid] = { sampler(i,j, sampleid, 0) , sampler(i,j, sampleid, 1) };

    double mse = calculate_mse_Gaussian( samples );
    result[ i + j*size] = mse;
    totalmse += mse;
  }

  //Writing the file
  //
  int res = stbi_write_hdr(filename.c_str(), size,size,1, result.data());

  std::cout<<"Mean mse: "<<totalmse/(double)(size*size)<<std::endl;
  return 0;
}
