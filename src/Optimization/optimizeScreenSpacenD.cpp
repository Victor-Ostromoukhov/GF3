// optimizeScreenSpace2D.cpp -- vo Oct 2020
// based on solidangle2D.cpp (May 2018, Oct 2019)
// conventions: we call iptsXY 2D integer points (t_int_point2D) in the screen space (pixels)
//              we call ptsUV 2D floating points (t_point2D) in the sampling space (a priory, U and V are independent)
//				we call iseedsUV integer seeds (t_int_point2D) needed to generate ptsUV through our OwenPlus
// action:  1. read input_tile_seeds_128x128.dat and feel iseedsUV[i] of size iseedsUV containing seeds for dimensions dim1 and dim2
//          2. genrate ptsUV[i], using our OwenPlus, with owen_tree_depth=32
//          3. apply simulated annealing: at each step, swap one pair of seeds, and test whether the cost function decreases
//          4. write newly optimized seeds into output_tile_seeds_128x128.dat (we write each 10000 steps; expected N of iterations is above 10.000.000
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <omp.h>
#include <assert.h>

#include <CLI11.hpp>
#include <Timer.hpp>
#include "Samplers/OwenScrambling.h"
#include "Samplers/SobolGenerator1D.h"
#include <SeedTileHelper.hpp>
#include "../Integration/Integration.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"
#include "../Samplers/SobolCascade.h"
#include "../Samplers/BlueTile.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#define TRACE_ITER 5000
#define EXPORT_ITER 5000
#define EXPORT_HDR_ITER 50000

typedef std::vector<double> Descriptor;
typedef unsigned int Index;

//Just a trace lambda
auto traceMultithread=[](uint *cpt, uint size){
  int tid = omp_get_thread_num();
  if (tid==0)
  {
    std::cout<<"\r"<<(*cpt)*100 << "               / "<<size*size <<std::flush;
    (*cpt)++;
  }
};

#define SEED 12345

//Preprocessing
Descriptor generateDescriptor(const std::vector<SobolGenerator1D<uint32_t> > &sobols,
	const SeedTile &seeds,
	const uint i,
	const uint j,
	const uint spp,
	const uint8_t Dimension) {
	std::vector< VecXDynamic > points;
	const int nbits = std::log2(spp);
	VecXDynamic sample(Dimension);
	std::vector<uint32_t> seed_array(Dimension);
	for (auto s = 0; s < spp; ++s) {
		// copie les seeds
		for (uint k = 0; k < Dimension; k++)
			seed_array[k] = seeds(i, j, k);
		// genere le sample
		getUniformND(sample.data(), sobols, seed_array.data(), Dimension, s, nbits, 32, true);
		points.push_back(sample);
	}
	Descriptor desc;
	for (auto k = 0; k < 16384; k++) {
		double v = calculate_mse(points, 2, 1, SEED + k);
		desc.push_back(v); //Gauss
	}
	return desc;
}// E_ij
// E_ij[k'] = 1 gausienne   != E_kl[k']
// TODO: VERIFIER SEED CALCULATE MSE



//Preprocessing
// l_2 distance between descriptors
double sqdirstDescriptor(const std::vector<Descriptor> &Descs,
                         const Index idx,
                         const Index idx2)
{
  double sum=0.0;
  double d=0.0;
  if(idx==idx2) return 0.0;
  for(auto i=0; i < Descs[idx].size(); ++i)
  {
    d =Descs[idx][i] - Descs[idx2][i] ;
    sum +=d*d;
  }
  return sum;
}



//Preprocessing
// pixel x pixel  distance between descriptors
void generateDescDistance(const std::vector<SobolGenerator1D<uint32_t> > &sobols,
                          const SeedTile &seeds,
                          std::vector<double>& descDistance,
                          std::vector<double>& mse,
                          uint size,
                          uint Dimension,
                          uint spp)
{
  long int s2= size*size;
  //DescriptorTile
  std::cout<<"Generating descriptors..."<<std::endl;
  std::vector<Descriptor> descriptorTile(size * size);
  uint cpt=0u;
#pragma omp parallel for
  for(auto idx=0; idx < size*size; ++idx)
  {
    if (idx%100 ==0) traceMultithread(&cpt,size);
    descriptorTile[idx] = generateDescriptor(sobols, seeds, idx % size, idx / size, spp, Dimension);
  }
  std::cout<<std::endl;
  std::cout<<"Ok"<<std::endl;
  
  std::cout<<"Computing mse per pixel"<<std::endl;
#pragma omp parallel for schedule(dynamic)
  for(auto idx=0; idx < s2; ++idx)
  {
    double sum=0.0;
    for(auto i=0; i < descriptorTile[idx].size(); ++i)
    {
      sum +=descriptorTile[idx][i]*descriptorTile[idx][i];
    }
    mse[idx] = sum / (double)descriptorTile[idx].size();
  }
  
  std::cout<<"Computing l2 distances between descriptors ("<< s2*s2 / 1024 / 1024<<"M)..."<<std::endl;
  // descDistance [ idx + idx2*size*szie ] contains the square
  // of the l2 distance between the descriptors of pixels idx and idx2
  
  cpt=0;
#pragma omp parallel for schedule(dynamic)
  for(long idx=0; idx < s2; idx++)
  {
    if (idx%100 ==0) traceMultithread(&cpt,size);
    for(auto idx2=idx; idx2 < s2; idx2++)
    {
      descDistance[idx  + idx2 * s2] = sqdirstDescriptor(descriptorTile, idx,idx2);
      descDistance[idx2 + idx * s2]  = descDistance[idx + idx2 * s2];
    }
  }
  std::cout<<std::endl;
  std::cout<<"Ok"<<std::endl;
}

//Preprocessing
// pixel x pixel  exp^(-distances(pix,pix)/ sigmap^2)
//(toroidal)
void generateSampleDistance(std::vector<double> &distances, uint size, double sigmaFactor, double sigma)
{
  double sigm2=sigma*sigma;
  long int s2=size*size;
  //Toroidal distance +
  auto sqdistsigma=[&](Index idx, Index idx2){ double dx = std::abs((double)(idx%size)-(idx2%size)), dy =std::abs((double)(idx/size)-(idx2/size));
    if (dx > (double)size/2.0)
      dx = (double)size - dx;
    if (dy >  (double)size/2.0)
      dy = (double)size - dy;
    return (dx*dx+dy*dy)/sigm2;};
  
  std::cout<<"Computing l2 distances between samples ("<< size * size * size * size / 1024 / 1024<<"M)..."<<std::endl;
  auto cpt=0u;
#pragma omp parallel for schedule(dynamic)
  for(long idx=0; idx < s2; idx++)
  {
    if (idx%100 ==0) traceMultithread(&cpt,size);
    for(auto idx2=idx; idx2 < s2; idx2++)
    {
      distances[idx  + idx2 * s2] = std::exp( - sqdistsigma(idx,idx2) );
      distances[idx2 + idx * s2]  = distances[idx + idx2 * s2];
    }
  }
  std::cout<<std::endl;
  std::cout<<"Ok"<<std::endl;
}

double energy(const std::vector<double> &sampleDistance,
              const std::vector<double> &DescDistance,
              const std::vector<Index> &currentPermutation,
              uint size)
{
  double sum=0.0;
  long int s2=size*size;
#pragma omp parallel for reduction(+:sum) schedule(dynamic)
  for(auto idx2=0; idx2 < s2; ++idx2)
    for(auto idx=0; idx < s2; ++idx)
  {
    assert(currentPermutation[idx] <s2);
    assert(currentPermutation[idx2] < s2);
    assert(currentPermutation[idx]  + currentPermutation[idx2]*s2 < s2*s2);
    sum += sampleDistance[idx + idx2*s2]  * DescDistance[ currentPermutation[idx]  + currentPermutation[idx2]*s2];
  }
  return sum;
}

double energyClamp(const std::vector<double> &sampleDistance,
                  const std::vector<double> &DescDistance,
                  const std::vector<Index> &currentPermutation,
                  const double sigmaFactor,
                   const double sigma,
                  uint size)
{
  double sum=0.0;
  int64_t s2=size*size;
#pragma omp parallel for reduction(+:sum) schedule(dynamic)
  for(auto idx=0; idx < s2; ++idx)
  {
    auto i = idx%size;
    auto j = idx/size;
    auto mini = (int)std::floor((double) i-sigmaFactor*sigma);
    auto maxi = (int)std::ceil((double)i+sigmaFactor*sigma);
    auto minj =  (int)std::floor((double)j-sigmaFactor*sigma);
    auto maxj = (int)std::ceil((double)j+sigmaFactor*sigma);
    for(auto jj=minj; jj<=maxj; ++jj)
    for(auto ii=mini; ii <= maxi ; ++ii)
    {
      auto idx2=((ii+size)%size)   +((jj+size)%size)*size;
      sum += sampleDistance[idx2 + idx*s2]  * DescDistance[ currentPermutation[idx2]  + currentPermutation[idx]*s2];
	}
  }
  return sum;
}


double energyClampLocal(const std::vector<double> &sampleDistance,
	const std::vector<double> &DescDistance,
	const std::vector<Index> &currentPermutation,
	const std::vector<Index> &perms,
	const double sigmaFactor,
	const double sigma,
	uint size)
{
	double sum = 0.0;
	int64_t s2 = size*size;
#pragma omp parallel for reduction(+:sum) schedule(dynamic)
	for (int64_t id = 0; id< perms.size(); id++)
	{
		auto idx = perms[id];
		auto i = idx%size;
		auto j = idx / size;
		auto mini = (int)std::floor((double)i - sigmaFactor*sigma);
		auto maxi = (int)std::ceil((double)i + sigmaFactor*sigma);
		auto minj = (int)std::floor((double)j - sigmaFactor*sigma);
		auto maxj = (int)std::ceil((double)j + sigmaFactor*sigma);
		for (auto jj = minj; jj<=maxj; ++jj)
			for (auto ii = mini; ii <= maxi; ++ii)
			{
				auto idx2 = ((ii + size) % size) + ((jj + size) % size)*size;
				sum += sampleDistance[idx2 + idx*s2] * DescDistance[currentPermutation[idx2] + currentPermutation[idx] * s2];
			}
	}
	for (int64_t id = 0; id < perms.size(); id++)
	{
		auto idx = perms[id];
		auto i = idx%size;
		auto j = idx / size;
		auto mini = (int)std::floor((double)i - sigmaFactor*sigma); 
		auto maxi = (int)std::ceil((double)i + sigmaFactor*sigma);
		auto minj = (int)std::floor((double)j - sigmaFactor*sigma);
		auto maxj = (int)std::ceil((double)j + sigmaFactor*sigma);
		if (maxi >= size) { maxi -= size; mini -= size; }
		if (maxj >= size) { maxj -= size; minj -= size; }
		for (int64_t id2 = id+1; id2 < perms.size(); id2++)
		{
			auto idx2 = perms[id2];
			auto i = idx2%size;		
			auto di = i - mini; if (di >= size ) di -= size ;
			if (di > maxi - mini) continue;

			auto j = idx2 / size;
			auto dj = j - minj; if (dj >= size ) dj -= size ;
			if (dj > maxj - minj) continue;
			sum -= sampleDistance[idx2 + idx*s2] * DescDistance[currentPermutation[idx2] + currentPermutation[idx] * s2];
		}
	}
	return sum;
}

void saveDescDistance(const std::vector<double> &descDistances,
                        const std::string &filename,
                        const size_t size)
{
  std::ofstream outfile(filename, std::ios::out);
  for(auto idx = 0; idx < size*size*size*size; ++idx)
    outfile << descDistances[idx]<<std::endl;
  outfile.close();
}

void loadDescDistance(std::vector<double> &descDistances,
                      const std::string filename,
                      const size_t size)
{
  double val;
  std::ifstream ifs(filename);
  for(auto idx=0; idx < size*size*size*size; ++idx)
  {
    ifs >> val;
    descDistances[idx] = val;
  }
  ifs.close();
}


void saveMSE(const std::vector<double> &mse,
             const std::string &filename,
                      const size_t size)
{
  std::ofstream outfile(filename, std::ios::out);
  for(auto idx = 0; idx < size*size; ++idx)
    outfile << mse[idx]<<std::endl;
  outfile.close();
}

void loadMSE(std::vector<double> &mse,
                      const std::string filename,
                      const size_t size)
{
  double val;
  std::ifstream ifs(filename);
  for(auto idx=0; idx < size*size; ++idx)
  {
    ifs >> val;
    mse[idx] = val;
  }
  ifs.close();
}



int main(int argc, char **argv)
{
  std::string in_fname, out_fname,dir_vectors_fname;
  uint8_t NbThreadsMax = omp_get_max_threads();
  uint Dimension=2;
  uint maxIter = 10000;
  uint seed= 1234556;
  uint batchSize=10;
  double kTemperature = 10.0;
  double sigma=2.1;
  double sigmaFactor = 2.5;
  uint spp=1;
  bool outputEnergy=false;
  uint iter=0;
  std::string exportDescDistance_fname="";
  std::string importDescDistance_fname="";
  std::string exportMSE_fname="";
  std::string importMSE_fname="";
  bool exitAfterExporting=false;
  std::string convertOutFile = "";
  std::vector<std::string> convertInFiles;
  CLI::App app { "optimizeScreenSpace2D" };
  app.add_option("-i", in_fname, "in_fname for tile_seeds, default: " + in_fname);// ->required()->check(CLI::ExistingFile);
  app.add_option("-o", out_fname, "out_fname for tile_seeds, default: " + out_fname );
  app.add_flag("--output-energy",outputEnergy, "export the energy to out-energy.dat" );
  app.add_option("-s,--spp", spp, "Number of samples default: " + std::to_string(spp) );
  app.add_option("-b,--batch", batchSize, "Batch size (number of swaps per iteration)");
  app.add_option("-M,--max", maxIter, "Max Iter" );
  app.add_option("--warmIter", iter, "Warm start (setting the iter)" );
  app.add_option("-k", kTemperature, "temperature decreasing factor, default: " + std::to_string(kTemperature) );
  app.add_option("--sigmaFactor", sigmaFactor, "Sigma factor to clamp the Gaussian, default: " + std::to_string(sigmaFactor) );
  app.add_option("--sigma", sigma, "Sigma for the spatial filter, default: " + std::to_string(sigma) );
  app.add_option("-n,--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + std::to_string(NbThreadsMax) );
  app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );//->required()->check(CLI::ExistingFile);;
  app.add_option("--dim", Dimension, "Dimension of the optimization ");
  app.add_option("--convertOutFile", convertOutFile, "Convert old tile to new tile ");
  app.add_option("--convertInFiles", convertInFiles, "Convert old tile to new tile ");
  
  app.add_option("--exportDescDistances", exportDescDistance_fname,"Compute the descriptor pairwise distance (and the MSE), export it to the given filename");
  app.add_option("--loadDescDistancesFromFile", importDescDistance_fname,"Import descriptor pairwise distance from file");
  app.add_option("--exportMSE", exportMSE_fname,"export the per pixel MSE -double-");
  app.add_option("--loadMSEFromFile", importMSE_fname,"Import MSE Per pixel from file");
  app.add_flag("--exitAfterExporting",exitAfterExporting, "exit after exporting MSE and DescDistances"   );

  CLI11_PARSE(app, argc, argv)
  
  std::cout<<"=============================================="<<std::endl;
  std::cout<<"Computing pex pixel integration:"<<std::endl;
  std::cout<<app.config_to_str(true,true);
  std::cout<<"=============================================="<<std::endl;
  
  
  //Fixing thread number
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions
  

  if (convertOutFile != "") {
	  SeedTile seedmap2tmp;
	  seedmap2tmp.loadTile(convertInFiles[0]);
	  BlueTile seedmap(seedmap2tmp.size, seedmap2tmp.dimension, convertInFiles.size());

	  for (int i = 0; i < convertInFiles.size(); i++) {
		  SeedTile seedmap2;
		  seedmap2.loadTile(convertInFiles[i]);
		  // blue tile, nouveau format
		  
		  // copier les valeurs
		  for (int y = 0; y < seedmap2.size; y++)
			  for (int x = 0; x < seedmap2.size; x++)
				  for (int d = 0; d < seedmap2.dimension; d++)
					  seedmap(x, y, d, i) = seedmap2(x, y, d);
	  }
	  // enregistre la nouvelle tuile 
	  seedmap.write(convertOutFile);
	  exit(0);
  }


  //Loading data
  SeedTile seedMap;
  seedMap.loadTile(in_fname);
  auto size = seedMap.size;


  
  if (seedMap.dimension != Dimension)
  {
    std::cout<<"[Error] the requested dimension and the dimension of the map do not match"<<std::endl;
    exit(2);
  }
  
  // E_ij = [  1spp intégration ] sur  16k gausiennes
  // exp( -d (ij - kl)^2/sigma) * || E_ij - E _kl||^2
  
  auto signstr=[](double d){ return (d>0)? "+": "-";};
  auto boolstr=[](bool d){ return (d)? "True": "False";};

  std::vector<SobolGenerator1D<uint32_t> > sobols;  // array of sobol data per dim
  loadSobolsFromFile(dir_vectors_fname, sobols);    // read sobols from file

  std::vector<double> descDistance( size * size * size * size,0.0);
  std::vector<double> msePerPixel( size * size,0.0);
  
  if (importDescDistance_fname=="")
    generateDescDistance(sobols, seedMap, descDistance,msePerPixel, size, Dimension, spp);
  else
    loadDescDistance(descDistance, importDescDistance_fname, size);
  
  if (importMSE_fname!="")
    loadMSE(msePerPixel, importMSE_fname,size);

  if (exportDescDistance_fname!= "")
  {
    saveDescDistance(descDistance,exportDescDistance_fname,size);
    if (exportMSE_fname!="")
      saveMSE(msePerPixel, exportMSE_fname, size);
  }
  
  if (exitAfterExporting)
    exit(0);
  
  std::vector<double> sampleDistance( size * size * size * size,0.0);
  generateSampleDistance(sampleDistance, size, sigmaFactor,sigma);
  
  std::mt19937_64 gen(std::time(nullptr));
  std::uniform_int_distribution<uint32_t> unif(0, size*size -1);
  std::uniform_real_distribution<double> u(0, 1.0);

  //Identity permutation to start with
  std::vector<Index> currentPermutation(size*size);
  for (auto i = 0; i < size*size; ++i) {
	  currentPermutation[i] = i;
  }
  
  CLI::Timer t("Energy");
  double initialEnergy = energy(sampleDistance,descDistance, currentPermutation,size);
  std::cout<<t<<std::endl;
  CLI::Timer tt("EneryC");
  double initialEnergyC = energyClamp(sampleDistance,descDistance, currentPermutation, sigmaFactor,sigma,size);
  std::cout<<tt<<std::endl;

  std::cout<<"Initial Energy= "<<initialEnergy<<"    "<<initialEnergyC<<std::endl;
  double currentEnergy = initialEnergyC;
  std::vector<float> result(size*size);
  std::string filename="descriptors-orig.hdr";
  double maxMSEVal = 0;
  for (int i = 0; i < size*size; i++) {
	  maxMSEVal = std::max(maxMSEVal, msePerPixel[i]);
  }
  for(auto i=0; i < size; ++i)
  {
    for(auto j=0; j < size; ++j)
    {
		double val = msePerPixel[i + j * size] / maxMSEVal;
		std::cout << val << " ";
      result[i+j*size] = val;
    }
  }
  std::cout << std::endl;
  stbi_write_hdr(filename.c_str(), size,size,1, result.data());
  
  
  auto SwapIndices=[&](size_t idx, size_t idx2){
	std::swap(currentPermutation[idx], currentPermutation[idx2]);
  };
  
  std::unordered_set<Index> idxSet;
  std::unordered_set<Index> idx2Set;
  std::vector<Index> perms;

  uint failure=0;
  uint currentBatchSize;
  double delta=0.0;
  
  
  while (iter < maxIter)
  {
    //batch of disinct pairs
    currentBatchSize = (uint)std::max(1, (int)batchSize - (int)failure);
	perms.resize(2 * currentBatchSize);
    idxSet.clear();
	int idp = 0;
	while (idxSet.size() < currentBatchSize) {
		uint32_t val = unif(gen);
		if (idxSet.find(val) == idxSet.end()) {
			idxSet.insert(val);
			perms[idp] = val; idp++;
		}
	}
    	
    idx2Set.clear();
    while (idx2Set.size() < currentBatchSize) { 
		auto val = unif(gen);
		if (idxSet.find(val) == idxSet.end()) {
			if (idx2Set.find(val) == idx2Set.end()) {
				idx2Set.insert(val);
				perms[idp] = val; idp++;
			}
		}
    }
    assert(idxSet.size() == currentBatchSize);
    assert(idxSet.size() == idx2Set.size());
	
	double newEnergyRemoved = energyClampLocal(sampleDistance, descDistance, currentPermutation, perms, sigmaFactor, sigma, size);

    //Do the swap in the batch
    for(auto it = idxSet.begin(), it2=idx2Set.begin();
        it != idxSet.end() ;
        ++it,++it2)
    {
      assert(*it != *it2);
      SwapIndices(*it,*it2);
    }
    //Check the new energy
	double newEnergyAdded = energyClampLocal(sampleDistance, descDistance, currentPermutation, perms, sigmaFactor, sigma, size);
	double newEnergy = currentEnergy - newEnergyRemoved*2 + newEnergyAdded*2;
    //double newEnergy = energyClamp(sampleDistance, descDistance, currentPermutation,sigmaFactor,sigma,size);
    double deltaE = (newEnergy - currentEnergy) / newEnergy;

    
    //Simulated Annealing
    bool accept = ( newEnergy > currentEnergy ) || ( u(gen)  > std::exp( - kTemperature / ((double)iter-1000.0*failure))); //deltaE / temperature));
    
    if ( accept )
    {
      currentEnergy = newEnergy;
      failure=0;
    }
    else
    {
      delta=0;
      
      //Do the swap back
      for(auto it = idxSet.begin(), it2=idx2Set.begin();
          it != idxSet.end() ;
          ++it,++it2)
      SwapIndices(*it,*it2);
      
      failure++;
    }
    
    
    if (iter % TRACE_ITER == 0)
    {
      std::cout<<iter<<" "<<currentEnergy<<" batchSize="
               <<currentBatchSize<<" failure="<<failure <<" "<< signstr(deltaE) <<" "<<boolstr(accept)<<" "
               <<" deltaE="<<deltaE <<" exp= "<< std::exp( - kTemperature / ((double)iter-1000.0*failure)) <<std::endl;
    
      if (outputEnergy)
      {
        std::ofstream outfile;
        outfile.open("out-energy.dat", std::ios_base::app); // append instead of overwrite
        outfile << iter<<" "<< currentEnergy<<std::endl;
        outfile.close();
      }
    }
    if (iter % EXPORT_ITER == 0)
    {
      std::cout<<"exporting "<<currentEnergy<<std::endl;
      SeedTile output( seedMap.size, seedMap.dimension);
      for(auto i=0; i < size; ++i)
      for(auto j=0; j < size; ++j)
      {
        Index idx = currentPermutation[ i + j*size];
        for(auto d=0; d < Dimension; ++d)
          output(i,j,d) = seedMap(idx % size, idx / size, d);
      }
      output.saveTile(out_fname);
      
      if (iter % EXPORT_HDR_ITER == 0)
      {
        std::vector<float> result(size*size);
        std::string filename="descriptors-"+std::to_string(iter/EXPORT_HDR_ITER)+".hdr";
        for(auto i=0; i < size; ++i)
        {
          for(auto j=0; j < size; ++j)
            result[i+j*size] = msePerPixel[ currentPermutation[ i + j*size]  ] / maxMSEVal;
        }
        stbi_write_hdr(filename.c_str(), size,size,1, result.data());
      }
      
    }
    iter++;
  }
  return 0;
}

