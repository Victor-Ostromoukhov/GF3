// owen.cpp
// vo 2004-2021

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include "CLI11.hpp"

//#include "OwenScrambling.h"
//#include "SobolGenerator1D.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)

//--------------------------------- structures
typedef uint32_t sobolInt;

//--------------------------------- global variables
double MSElimit = 3.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];

int main(int argc, char **argv) {
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i){
        powers_of_two[i] = 1u << i;
    }
	CLI::App app { "owen" };
	unsigned int nDims=4, ntrials = 1, firstDim = 0, npts = 16;
	std::string filename = "out.dat";
	std::string dir_vectors_fname = "data/sobol_init_tab.dat";
	bool dbg_flag = false;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    uint32_t imposed_owen_tree_depth = 0;
    uint32_t starting_from = 0;

	uint32_t seed = time(NULL);
	app.add_option("-n,--npts", npts, "number of points to generate, default: "+ std::to_string(npts)+")" );
	app.add_option("--nd,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims)+")" );
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename + ")" );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: 13374269");
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );
	app.add_option("-f,--firstDim", firstDim, "d dimensions, starting from firstDim, default: "+ std::to_string(firstDim)+")" );
	app.add_option("-p,--permut", owen_permut_flag, "permut or not? when owen_permut_flag=0 ==> Sobol, default: "+ std::to_string(owen_permut_flag)+")" );
	app.add_option("-m,--max_tree_depth_32_flag", max_tree_depth_32_flag, "Owen or OwenPlus ? default: OwenPlus max_tree_depth_32_flag="+ std::to_string(max_tree_depth_32_flag)+"" );
	app.add_option("-t,--tree_depth", imposed_owen_tree_depth, "real tree dephth in Owen+" );
	app.add_option("--start", starting_from, "N-D points, starting from this count, default: "+ std::to_string(starting_from)+")" );
	app.add_option("--ntrials", ntrials, "number of pointsets to generate, default: "+ std::to_string(ntrials)+")" );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)+")" );
	CLI11_PARSE(app, argc, argv)

	std::ofstream fout(filename, std::ofstream::out);
	std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    loadSobolsFromFile(dir_vectors_fname, sobols, seed, dbg_flag);		// read sobols from file
    if(dbg_flag) {
	    for(int i = 0; i < nDims; i++) {
	        cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
	        for(int j = 0; j < sobols[i].s; j++) cout << " " << sobols[i].m[j];
	        cout  << endl;
	    }
	    cout  << " after filling 32 levels: " << endl;
	    for(int i = 0; i < nDims; i++) {
	        cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
	        for(int j = 0; j < 32; j++) cout << " " << sobols[i].m[j];
	        cout  << endl;
	    }
    }

	uint32_t owen_tree_depth = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
	if(imposed_owen_tree_depth > 0) owen_tree_depth = imposed_owen_tree_depth;
	if(max_tree_depth_32_flag) owen_tree_depth = OWENPLUS_TREE_DEPTH;

	for (int itrial = 0; itrial < ntrials; itrial++) {
    	uint32_t ipt_from = starting_from, ipt_to = starting_from + npts, ipt_step = 1;
		for (int ipt = ipt_from; ipt < ipt_to; ipt += ipt_step) {
			if (dbg_flag) cout << "pt = " << ipt << endl;
			for (unsigned int idim = firstDim; idim < firstDim+nDims; idim++) {
				float res_float = getOwenPlus1D(sobols, idim, ipt, owen_tree_depth, dbg_flag, owen_permut_flag);
				fout << std::setprecision(20) << std::fixed << std::setprecision(20) << res_float << " ";
			}
			fout << endl;
		}
		if (itrial != ntrials-1) fout << "#" << endl;
	}

	fout.close();

	if (dbg_flag) cout << nDims <<"-d points written into " << filename << " : " << npts*ntrials << " pts." << endl;

}
