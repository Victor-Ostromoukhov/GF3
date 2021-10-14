
#ifndef _SEEDTILE_H
#define _SEEDTILE_H

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <cassert>
#include <cstdint>

//Class to define and handle tile of nD seeds.
struct SeedTile{
  
  SeedTile() = default;
  
  
  /// Constructor
  /// @param aSize size of the tile in screenspace
  /// @param aDimension the dimension of the seedmap (e.g. 2D = UV)
  SeedTile(const uint32_t aSize, const uint32_t aDimension)
  {
    size = aSize;
    dimension = aDimension;
    this->seedmap.resize(size*size*dimension);
  }
  
  uint32_t &operator()(const uint32_t i, const uint32_t j, const uint32_t dim) { assert(dim < dimension); assert((i+j*size) + dim*size*size < seedmap.size()); return seedmap[(i+j*size) + dim*size*size];}
  const uint32_t &operator()(const uint32_t i, const uint32_t j, const uint32_t dim) const { assert(dim < dimension); assert((i+j*size) + dim*size*size < seedmap.size()); return seedmap[i+j*size + dim*size*size];}
 
  uint32_t &operator()(const uint32_t idx) {return seedmap[idx];}
  const uint32_t &operator()(const uint32_t idx) const {return seedmap[idx];}
  
  void saveTile(const std::string &filename)
  {
    std::ofstream ofs(filename,  std::ios::out);
    ofs << size<< " "<<dimension <<std::endl;
    for(auto idx=0; idx < size*size*dimension; ++idx)
       ofs << this->operator()(idx)<<" ";
    ofs <<std::endl;
    ofs.close();
  }
  
  void loadTile(const std::string &filename)
  {
    std::ifstream ifs(filename);
    assert(ifs.is_open());
    uint32_t val;
    ifs >> size;
    ifs >> dimension;
    seedmap.resize(size*size*dimension);
    std::cerr<<"Loading tilemap " << filename << " (size="<<size<<", dimension="<<dimension<<")"<<std::endl;
    for(auto idx=0; idx < size*size*dimension; ++idx)
    {
      ifs >> val;
      this->operator()(idx) = val;
    }
    ifs.close();
  }
  
  uint32_t dimension; ///Dimension of the Seed map
  uint32_t size;   ///Size of the tile
  std::vector<uint32_t> seedmap;  ///The data
};

#endif

