#include <iostream>

#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp.hpp"
#include "samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp.hpp"


int main()
{
  std::cout<<"The first 3d sample in pixel (0,0):  "<<std::endl;
  for(auto sampleid=0; sampleid<1; sampleid++)
  {
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_1spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    
    for(auto i=0;i < 3; ++i)
       std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_2spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
      std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_4spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_8spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_16spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_32spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_64spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_128spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
    for(auto i=0;i < 3; ++i)
    std::cout<< samplerBlueNoiseErrorDistribution_128x128_OptimizedFor_2d2d2d2d_256spp(0, 0, sampleid, i)<<" ";
    std::cout<<std::endl;
  }
  return 0;
}
