#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <algorithm>    // std::min
#include <math.h>       /* pow */
#include "CLI11.hpp"

#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"
#include "../Tools/io.h"


using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)

//--------------------------------- structures
typedef uint32_t sobolInt;
typedef struct t_int_point2D {int	x ; int	y ;} t_int_point2D;
typedef struct t_point2D {double	x ; double	y ;} t_point2D;

//--------------------------------- global variables
double MSElimit = 3.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];

int main(int argc, char **argv) {
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i){
        powers_of_two[i] = 1u << i;
    }
	unsigned int nDims=2, ntrials = 1000000000, owen_tree_depth=32;
	std::string filename = "out.dat";
	std::string dir_vectors_fname = "data/vo_sobol_init_tab.dat";
	bool dbg_flag = false;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    int nSectors = 16, octaveMax = 6, dim1 = 3, dim2 = -1;
    double scalePow = 1.5;
	uint32_t seed = time(NULL);
	uint8_t NbThreadsMax = omp_get_max_threads();

	CLI::App app { "findGoodScreenSpaceSeeds" };
	app.add_option("--nd,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims)+")" );
	app.add_option("-i,--dim1", dim1, "dim1, default: "+ std::to_string(dim1)+")" );
	app.add_option("-j,--dim2", dim2, "dim2, default: "+ std::to_string(dim1+1)+")" );
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename + ")" );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: 13374269");
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );
	app.add_option("--ntrials", ntrials, "number of pointsets to generate, default: "+ std::to_string(ntrials)+")" );
	app.add_option("-m,--octaveMax", octaveMax, "max octave value, default: " + std::to_string(octaveMax) );
	app.add_option("-n,--nSectors", nSectors, "factor of subdivision of the domain [0,1], per dimension , default: " + std::to_string(nSectors) );
	app.add_option("-p,--pow", scalePow, "scalePow in scaleBary = pow((double)npts, scalePow); , default: " + std::to_string(nSectors) );
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + std::to_string(NbThreadsMax) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)+")" );
	CLI11_PARSE(app, argc, argv)
    if(dim2 == -1)  dim2 = dim1 + 1;
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);	// uniform distribution of integers between 0 and 2^32-1

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

//    // double-testing
//	std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
//    loadSobolsFromFile(dir_vectors_fname, sobols, seed, dbg_flag);		// read sobols from file
//    for(int i = 0; i < sobols[1].m.size(); i++ )
//    	cout << std::bitset<32>(sobols[1].m[i]) << endl;
//    for(int i = 0; i < sobols[2].m.size(); i++ )
//    	cout << std::bitset<32>(sobols[2].m[i]) << endl;
//	for (int ipt = 0; ipt < 16; ipt++) {
//		t_point2D pts;
//		pts.x = getOwenPlus1D_with_seed(sobols, dim1, ipt, 4244163743 );	// 4244163743 1117257301 -> 0.3179032207 0.5295686722
//		pts.y = getOwenPlus1D_with_seed(sobols, dim2, ipt, 1117257301 );	// 3025097427 4293809212 -> 0.9502353668 0.5428119302
//		cout << 16 * getOwenPlus1D(sobols, 1, ipt, owen_tree_depth, dbg_flag, false) << " \t" <<
//				16 * getOwenPlus1D(sobols, 2, ipt, owen_tree_depth, dbg_flag, false) << " -> \t" << pts.x << " " << pts.y << endl;
//	}
//	exit(0);

	std::vector<SobolGenerator1D<uint32_t> > init_sobols;	// array of sobol data per dim
    loadSobolsFromFile(dir_vectors_fname, init_sobols, seed, dbg_flag);		// read sobols from file

#pragma omp parallel
    while (true) {
    	std::vector<SobolGenerator1D<uint32_t> > sobols = init_sobols;
        uint32_t seedX = unif32bits(gen);
        uint32_t seedY = unif32bits(gen);
        if(seedX == seedY) continue;
       	vector<t_int_point2D> sectors, scaleddBarycenters;
    	sectors.clear();
    	scaleddBarycenters.clear();
    	t_point2D barycenter, pt, scaleddBary, ref_barycenter;
        for (int ioctave = 0; ioctave <= octaveMax ; ioctave += 2) {
        	int npts = powers_of_two[ioctave];
        	double scaleBary = pow((double)npts, scalePow);
        	t_int_point2D sector;
         	barycenter.x = barycenter.y = 0.;
    		for (int ipt = 0; ipt < npts; ipt++) {
    			pt.x = getOwenPlus1D_with_seed(sobols, dim1, ipt, seedX, 32, dbg_flag, true);
    			pt.y = getOwenPlus1D_with_seed(sobols, dim2, ipt, seedY, 32, dbg_flag, true);
//    			cout << ioctave << " " << ipt << " | " << seedX << " " << seedY << " -> " << pt.x << " " << pt.y << endl;
    			barycenter.x += pt.x;
    			barycenter.y += pt.y;
    		}
			barycenter.x /= (double)npts;
			barycenter.y /= (double)npts;
			scaleddBary.x = scaleBary*(barycenter.x - .5) + .5;
			scaleddBary.y = scaleBary*(barycenter.y - .5) + .5;
			sector.x = round((double)nSectors * scaleddBary.x);
			sector.y = round((double)nSectors * scaleddBary.y);
			if(ioctave == 0) ref_barycenter = barycenter;
			// clump it
			sector.x = std::min(std::max(sector.x,0),nSectors);
			sector.y = std::min(std::max(sector.y,0),nSectors);
			sectors.push_back(sector);
			scaleddBarycenters.push_back(sector);
        }
        bool found = true;
        for (int ioctave = 0; ioctave < octaveMax/2 ; ioctave += 1) {
        	found = found && ((sectors[ioctave].x == sectors[ioctave+1].x) && (sectors[ioctave].y == sectors[ioctave+1].y ));
        }

        if(found)
#pragma omp critical (output)
        	{
            	cout << sectors[0].x << " " << sectors[0].y << " " << ref_barycenter.x << " " << ref_barycenter.y << " " << seedX << " " << seedY << " " << endl;
            	for (int ioctave = 0; ioctave <= octaveMax/2 ; ioctave += 1) {
            		cout << ioctave << " " << scaleddBarycenters[ioctave].x << " " << scaleddBarycenters[ioctave].y << endl;
            	}
            	cout.flush();
        	}
    }
}
