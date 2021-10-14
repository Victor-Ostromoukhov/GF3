// integrateSobolPlusPlus2D.cpp

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <math.h>       /* fmod */
#include "CLI11.hpp"

#include "../Integration/Integration.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"
#include "../Tools/io.h"
#ifdef _MSC_VER 
#include "../Tools/drand48.h"
#endif

using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)

//--------------------------------- structures
typedef uint32_t sobolInt;

//--------------------------------- global variables
double MSElimit = 3.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];


//--------------------------------- routines related to integration
//template<int DIM>
void integratePointSets(
        const std::vector<SobolGenerator1D<sobolInt> >& sobols,
        const unsigned int firstDim,
        const unsigned int nDims,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        const unsigned int octaveStep,
        const unsigned int nTrials_integration,
		mt19937_64& gen,
        std::ofstream& out,
		unsigned int shift_switch,
        const bool dbg_flag,
        const int integrandType = 1,
		bool owen_permut_flag = true,
		bool max_tree_depth_32_flag = true
		) {

    std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);	// uniform distribution of integers between 0 and 2^32-1
    std::uniform_real_distribution<double> unifreal01(0, 1.0);	// uniform distribution of reals between 0 and 1.

    vector<double> integral_estimations;
    vector<t_pointND> mse_data;
    double mse_accumulator;
	sobolInt IDcode;
    uint32_t owen_seeds[32];	// 32 dims max
    uint32_t digital_shift1, digital_shift2;
    double CP_shift1, CP_shift2;

	for (int ioctave = octaveMin; ioctave <= octaveMax; ioctave += octaveStep) {

		size_t npts = powers_of_two[ioctave];
        vector<VecX<2>> points2(npts);
        vector<VecX<3>> points3(npts);
        vector<VecX<4>> points4(npts);
        vector<VecX<5>> points5(npts);
        vector<VecX<6>> points6(npts);
        vector<VecX<8>> points8(npts);
        vector<VecX<10>> points10(npts);
        vector<VecX<12>> points12(npts);
    	mse_accumulator = 0.;
//#pragma omp parallel for reduction(+:mse_accumulator) // works poorly on ubu12, ubu32
		for (int iset = 0; iset < nTrials_integration; iset++) {

			switch(shift_switch) {
			    case 1 :
			    	for (int idim = 0; idim < nDims; idim++)
			            owen_seeds[idim] = unif32bits(gen);
	                break;
			    case 2 :
		            digital_shift1 = unif32bits(gen);
		            digital_shift2 = unif32bits(gen);
	                break;
			    case 3 :
		            CP_shift1 = unifreal01(gen);
		            CP_shift2 = unifreal01(gen);
	                break;
			}
			for (sobolInt ipt = 0; ipt < npts; ipt++) {
				switch(shift_switch) {
				    case 1 :
				    	for (int idim = 0; idim < nDims; idim++) {
							IDcode = sobols[firstDim+idim].getSobolInt(ipt);
							if(owen_permut_flag) {
								uint32_t owen_tree_depth;
								if(max_tree_depth_32_flag)
									owen_tree_depth = OWENPLUS_TREE_DEPTH;
								else
									owen_tree_depth = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
								IDcode = OwenScrambling(IDcode, owen_seeds[idim], owen_tree_depth);
							}
							switch(nDims) {
				    		case 2:
				    			points2[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 3:
				    			points3[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 4:
				    			points4[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 5:
				    			points5[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 6:
				    			points6[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 8:
				    			points8[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 10:
				    			points10[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
				    		case 12:
				    			points12[ipt][idim] = ((double) IDcode / (double) UINT32SOBOLNORM);
				    			break;
							}
				    	}
		                break;
				    case 2 :
						IDcode = sobols[firstDim+0].getSobolInt(ipt);
						IDcode ^= digital_shift1;
						points2[ipt][0] = ((double) IDcode / (double) UINT32SOBOLNORM);
		                IDcode = sobols[firstDim+1].getSobolInt(ipt);
						IDcode ^= digital_shift2;
						points2[ipt][1] = ((double) IDcode / (double) UINT32SOBOLNORM);
		                break;
				    case 3 :
						IDcode = sobols[firstDim+0].getSobolInt(ipt);
						double x = ((double) IDcode / (double) UINT32SOBOLNORM) + CP_shift1;
						points2[ipt][0] = std::fmod(x, 1.);
		                IDcode = sobols[firstDim+1].getSobolInt(ipt);
		                double y = ((double) IDcode / (double) UINT32SOBOLNORM) + CP_shift2;
		                points2[ipt][1] = std::fmod(y, 1.);
		                break;
				}


			}
			double this_mse;
			switch(nDims) {
    		case 2:
    			this_mse = calculate_mse(points2, integrandType);
    			break;
    		case 3:
    			this_mse = calculate_mse(points3, integrandType);
    			break;
    		case 4:
    			this_mse = calculate_mse(points4, integrandType);
    			break;
    		case 5:
    			this_mse = calculate_mse(points5, integrandType);
    			break;
    		case 6:
    			this_mse = calculate_mse(points6, integrandType);
    			break;
    		case 8:
    			this_mse = calculate_mse(points8, integrandType);
    			break;
    		case 10:
    			this_mse = calculate_mse(points10, integrandType);
    			break;
    		case 12:
    			this_mse = calculate_mse(points12, integrandType);
    			break;
			}
			mse_accumulator += this_mse;
		}
		double mean_mse = mse_accumulator/(double)nTrials_integration;
		cout << "integrandType=" << integrandType << " =====================> " << firstDim << " nDims=" << nDims << " " << npts << " "  << nTrials_integration << " -> " << " " << mean_mse << endl;
		out  << std::setprecision(20) << npts << " \t" << mean_mse << endl;
	}
} // integratePointSets


int main(int argc, char **argv) {
	srand48( time(NULL) );
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i){
        powers_of_two[i] = 1u << i;
    }
	unsigned int firstDim=0, nDims=2;
	unsigned int octaveMin = 0;
	unsigned int octaveMax = 14;
	unsigned int nTrials_integration = 32;
	unsigned int shift_switch  = 1;	// 1 : Owen+  2 : digital shift 3 : CP

	std::string filename = "out.dat";
	std::string dir_vectors_fname = "data/sobol_init_tab.dat";
    uint8_t NbThreadsMax = omp_get_max_threads();

	bool dbg_flag = false;
    bool owenPlus_permut_flag = true;
    int integrandType = 2;
	uint32_t seed = 13374269;
	bool owen_permut_flag = true, max_tree_depth_32_flag = true;

	CLI::App app { "integrateSobolPlusPlus2D" };
	app.add_option("-n,--nDims", nDims, "number of dimensions default: " + std::to_string(nDims) );
	app.add_option("-f,--firstDim", firstDim, "Index of the first Sobol pair (==firstDim) convetion: firstDim=0 : van der Corput, default: " + std::to_string(firstDim) ); //->required();
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename + ")" );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: 13374269");
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );
	app.add_option("-p,--permut", owen_permut_flag, "permut or not? when owen_permut_flag=0 ==> Sobol, default: "+ std::to_string(owen_permut_flag)+")" );
	app.add_option("-m,--max_tree_depth_32_flag", max_tree_depth_32_flag, "Owen or OwenPlus ? default: OwenPlus max_tree_depth_32_flag="+ std::to_string(max_tree_depth_32_flag)+"" );
	app.add_option("--octaveMin", octaveMin, "min octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMin) );
	app.add_option("--octaveMax", octaveMax, "max octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMax) );
	app.add_option("--nTrials_integration", nTrials_integration, "nTrials_integration, default: " + std::to_string(nTrials_integration)+")");
	app.add_option("--owenPlus_permut_flag", owenPlus_permut_flag, "owenPlus_permut_flag, default: true");
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads (def: max number of thread of the OS = " + std::to_string(NbThreadsMax)+")");
	app.add_option("--shift", shift_switch, " type of shift -- 1 : Owen+  2 : digital shift 3 : CP  (default= " + std::to_string(shift_switch)+")");
	app.add_option("--integrandType", integrandType, "integrandType, possible values : 1 Heaviside 2 Gauss 3 Smooth 4 Cont 5 HeaviBell 6 HeaviCont 7 HeaviGauss default: " + std::to_string(integrandType) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: false");
	CLI11_PARSE(app, argc, argv)

	std::ofstream out(filename, std::ofstream::out);

    std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    if(dbg_flag) cout << "Loading Sobol dir_vectors from " << dir_vectors_fname << endl;
    loadSobolsFromFile(dir_vectors_fname, sobols, seed);		// read sobols from file


    std::uniform_int_distribution<uint32_t> unifFull(0);
    std::mt19937_64 gen(seed + firstDim);

    cout << dir_vectors_fname << " integrateSobolPlusPlusGauss2D integrandType=" << integrandType << " firstDim=" << firstDim << " nDims=" << nDims << " for npts between 2^1 and 2^" << octaveMax << " output into " << filename << endl;

    integratePointSets(sobols, firstDim, nDims, octaveMin, octaveMax, 1, nTrials_integration, gen, out, shift_switch, dbg_flag, integrandType, owen_permut_flag, max_tree_depth_32_flag);
    out.close();
}
