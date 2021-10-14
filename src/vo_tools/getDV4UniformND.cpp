// getDV4UniformND.cpp
// vo March 2021

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <ctime>
#include <omp.h>
#include "CLI11.hpp"

//#include "OwenScrambling.h"
//#include "SobolGenerator1D.h"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

#include "../Samplers/SobolCascade.h"
#include "../Samplers/BlueTile.h"

using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)
//#define PURE_SOBOLO_LIMIT_DIM	100		// below this limit: UniformND, above : pure Sobol/Owen

#define TILE_SIZE	1
#define NDIMENSIONS 200
#define NSPP 1

//--------------------------------- structures
typedef uint32_t sobolInt;


struct type_pointND {double * coord ; };
struct type_pointND_int {uint32_t * coord ; };

//--------------------------------- global variables
struct t_sobol {
    int d, s, a;
    std::array<uint32_t, CHAR_BIT*sizeof(uint32_t)> m;
//    std::vector<uint32_t> m;
};
struct t_stats {
    int Nbpts;
    float Mean, Var, Min, Max, Analytical, MSE, NbPtsets, k;
};

#include "sobol_init_tab_joe_kuo.h"

uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];

double discrepancy_refs[101][19] =
#include "discrepancy_refs_OwenPureFirstDim2.dat"
;
//--------------------------------- routines related to GeneralizedL2Discrepancy Formula from Lemieux p185
inline double calculate_GeneralizedL2Discrepancy(const unsigned int nDims, const unsigned int N, double **points )
{
	unsigned int D = nDims;
	long double a, factor_b, factor_c;
	a = pow((4.0/3.0), D);
	factor_b = 2.0/(double)N;
	factor_c = 1.0/(double)(N);
	factor_c *= factor_c;
	long double sumb = 0.0;
#pragma omp parallel for reduction(+:sumb)
	for(unsigned int i=0; i<N; i++)
	{
		long double prodb = 1.0;
		for(unsigned int j=0; j<D; j++)
		{
			double uij = points[i][j];
			prodb *= ((3.0 - uij*uij)/2.0);
		}
		sumb += prodb;
	}
	long double sumc = 0.0;
#pragma omp parallel for reduction(+:sumc)
	for(uint i=0; i<N; i++)
	for(unsigned int iprime=0; iprime<N; iprime++)
	{
		long double prodc = 1.0;
		for(unsigned int j=0; j<D; j++)
		{
			double uij = points[i][j];
			double uiprimej = points[iprime][j];
			double m = uij > uiprimej ? uij : uiprimej;//std::max(uij, uiprimej);
			prodc *= (2.0 - m);
		}
		sumc += prodc;
	}
	long double tmp0 = factor_b*sumb;
	long double tmp1 = factor_c*sumc;
	return sqrtl(a -tmp0 + tmp1);
} // calculate_GeneralizedL2Discrepancy

void writeStats(const std::string& out_discrepancy_fname,
		const vector<t_stats> stats,
	    const bool dbg_flag = false ) {
    std::ofstream out_file(out_discrepancy_fname);
    if (!out_file.is_open()) {
        std::cerr << out_discrepancy_fname << " cannot open." << std::endl;
        exit(1);
    };
	out_file << "#Nbpts	#Mean	#Var	#Min	#Max	#Analytical	#MSE	#NbPtsets	#Nbintegrands\n" ;
	for (int i = 0; i < stats.size(); i++) {
		out_file << stats[i].Nbpts << " " << stats[i].Mean << " " << stats[i].Var << " " << stats[i].Min << " "
			<< stats[i].Max << " " << stats[i].Analytical << " " << stats[i].MSE << " " << stats[i].NbPtsets << " " << stats[i].k << endl;
	}
    if (dbg_flag) std::cout << "writing into file " << out_discrepancy_fname << " " << stats.size() << " entries." << std::endl;
	out_file.close();
}

inline void writeSobols(
        const std::string& out_dir_vectors_fname,
		const std::vector<SobolGenerator1D<uint32_t> >& sobols,
	    const bool dbg_flag = false ) {
    std::ofstream tableFile(out_dir_vectors_fname);
    if (!tableFile.is_open()) {
        std::cerr << out_dir_vectors_fname << " cannot open." << std::endl;
        exit(1);
    };
	for(int i = 1; i < sobols.size(); i++) {	// we don't write entry 0 for van-der-corput
		tableFile << i // instead of sobols[i].d
				<< " " << sobols[i].s << " " << sobols[i].a << "\t" ;
		for(int j = 0; j < sobols[i].s; j++)
			tableFile << " " << sobols[i].m[j];
		tableFile << std::endl;
	}
    if (dbg_flag) std::cout << "writing into file " << out_dir_vectors_fname << " " << sobols.size() << " entries." << std::endl;

	tableFile.close();
}	// writeSobols

int main(int argc, char **argv) {
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i) powers_of_two[i] = 1u << i;
	unsigned int nDims = 6, ntrials = 4000000000, firstDim = 0, npts = 16;
	unsigned int tstDim = -1;
	std::string output_dv = "out_dv.dat";
	std::string output_discrepancy = "out_discrepancy.dat";
	std::string input_dir_vectors = "data/vo_sobol_init_tab.dat";
	bool dbg_flag = false;
    bool owen_permut_flag = true;
    uint32_t owen_tree_depth = 32, add_owen_tree_depth = 12;
    uint32_t nbits;
    uint32_t ocatve_from = 6, octave_to = 18, octave_step = 2;
    double kmax_limit = 1.0;
	uint8_t NbThreadsMax = omp_get_max_threads();
	double rand_val = 0.;
	uint32_t seed = time(NULL);

	CLI::App app { "owen" };
	app.add_option("-d,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims) );
	app.add_option("-t,--testDim", tstDim, "ind of the test dimension, default: "+ std::to_string(tstDim)+")" );
	app.add_option("--odv,--output_dv", output_dv, "output dir_vectors filename (ascii file), default: " + output_dv  );
	app.add_option("--odiscrepancy,--output_discrepancy", output_discrepancy, "output discrepancy filename (ascii file), default: " + output_discrepancy );
	app.add_option("--id,--idv", input_dir_vectors, "input dir_vectors (ascii file), default: " + input_dir_vectors );
    app.add_option("--of,--ocatve_from", ocatve_from, "starting octave, default: "+ std::to_string(ocatve_from)  );
    app.add_option("--ot,--ocatve_to", octave_to, "ending octave, default: "+ std::to_string(octave_to) );
    app.add_option("--os,--octave_step", octave_step, "octave step, default: " + std::to_string(octave_step)  );
	app.add_option("-p,--permut", owen_permut_flag, "permut or not? when owen_permut_flag=0 ==> Sobol, default: "+ std::to_string(owen_permut_flag) );
//	app.add_option("-a,--add_owen_tree_depth", add_owen_tree_depth, "additional tree_depth wrt std owens depth, default: "+ std::to_string(add_owen_tree_depth) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)  );
    app.add_option("--ntrials", ntrials, "number of pointsets to generate, default: "+ std::to_string(ntrials) );
    app.add_option("-k,--kmax", kmax_limit, "limit for kmax, default: " + std::to_string(kmax_limit) );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: " + std::to_string(seed) );
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + std::to_string(NbThreadsMax) );
	CLI11_PARSE(app, argc, argv)

	if(tstDim == -1) tstDim = nDims;
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

	//========================== allocation of tile_seeds
	std::mt19937_64 gen( seed );
    BlueTile seeds(TILE_SIZE, NDIMENSIONS, NSPP);
	//========================== end of allocation of tile_seeds

	std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    std::vector<uint32_t> d;
    std::vector<uint32_t> s;
    std::vector<uint32_t> a;
    std::vector<std::vector<uint32_t>> m;
    std::ifstream tableFile(input_dir_vectors);
    if (!tableFile.is_open()) {
        std::cerr << input_dir_vectors << " cannot be read." << std::endl;
        exit(1);
    };
    load_init_table(tableFile, d, s, a, m);
    std::cout << "loading file " << input_dir_vectors << " " << m.size() << " entries." << std::endl;
    sobols.resize(m.size());

    double **points;
    points = (double**)malloc(sizeof(double*) * powers_of_two[octave_to] );
	for (unsigned int i = 0; i < powers_of_two[octave_to]; i++)
		points[i] = (double *)malloc(sizeof(double) * nDims );
    int found = 0;
    vector<t_stats> stats;
	for (unsigned int itrial = 1; itrial <= ntrials; itrial++) {
	    double kmax = 0;
	    if(dbg_flag) {
		    for(int i = 0; i < nDims; i++) {
		        cout << d[i] << " " << s[i] << " " << a[i] << " ";
		        for(int j = 0; j < s[i]; j++) cout << " " << m[i][j];
		        cout  << endl;
		    }
	    }
        d[nDims-1] = d[nDims-1];	// <<<<<<   +0/+1/+2 now, it's not the dim index, but the length of effective dv
        s[nDims-1] = s[tstDim-1];	//sobolDataJoeKuo[tstDim-2].s;
        a[nDims-1] = a[tstDim-1];	//sobolDataJoeKuo[tstDim-2].a;
        m[nDims-1].resize(s[nDims-1]);
        m[nDims-1][0] = 1;
        for(int j = 1; j < s[nDims-1]; j++)  {
            std::uniform_int_distribution<uint32_t> unif32bits(0, powers_of_two[j-1]);	// uniform distribution of uint32_ts between 0 and 2^32-1
        	m[nDims-1][j] = 1 + 2 * unif32bits(gen);	// sobolDataJoeKuo[tstDim-2].m[j];
        }
	    init_sobols(sobols, d, s, a, m, seed, false);
	    if(dbg_flag) {
	        cout  << "after assignment of " << nDims-1 << "-th entry" << endl;
	    	int i = nDims-1;
	        cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
	        for(int j = 0; j < sobols[i].s; j++) cout << " " << sobols[i].m[j];
	        cout  << endl;
	    }
	    stats.clear();
		for (int octave = 2; octave <= octave_to; octave += octave_step) {
		    seeds.random_init(gen);         // init par generateurs c++
		    struct t_stats this_stats;
			npts = powers_of_two[octave];
			nbits = octave;
			for (int ipt = 0; ipt < npts; ipt++)
				getUniformND(points[ipt] , sobols, seeds(0, 0, 0), nDims, ipt, nbits, owen_tree_depth, owen_permut_flag);
			double k, discrepancy = calculate_GeneralizedL2Discrepancy(nDims, npts, points);
			k = discrepancy / discrepancy_refs[nDims][octave];
			time_t now = time(0);
			tm *ltm = localtime(&now);
			if (octave >= ocatve_from) kmax = max(k, kmax);
			if( ( ( octave == octave_to - octave_step) && (kmax < kmax_limit) ) || dbg_flag)
				cout << ltm->tm_hour<<":"<<ltm->tm_min<<":"<<ltm->tm_sec << " npts = " << npts << " \t d -> " << discrepancy  << " \tk -> " << k << " \t kmax -> " << kmax << " \t kmax_limit -> " << kmax_limit << endl;
			if (kmax > kmax_limit && octave == octave_to ) cout << ltm->tm_hour<<":"<<ltm->tm_min<<":"<<ltm->tm_sec << " npts = " << npts << "--------- nope. " << kmax << endl;
			if (kmax > kmax_limit && octave >= ocatve_from) {
				if(dbg_flag) cout << ltm->tm_hour<<":"<<ltm->tm_min<<":"<<ltm->tm_sec << " npts = " << npts << "--------- nope. " << kmax << endl;
				break;
			}
			this_stats.Nbpts = npts;
			this_stats.Mean = discrepancy;
			this_stats.Var	=	this_stats.Min	=	this_stats.Max	=	this_stats.Analytical	=	this_stats.MSE	=	0;
			this_stats.NbPtsets	= 1;
			this_stats.k = k;
			stats.push_back ( this_stats );
			if(octave == octave_to) {
				output_dv = "data_discrepancy/" + to_string(nDims) + "D/optim_dv/" + to_string(kmax) + "_" + to_string(tstDim) + "_";
				for(int j = 0; j < s[nDims-1]; j++)
					if(j == s[nDims-1] - 1)
						output_dv = output_dv + to_string(m[nDims-1][j]) ;
					else
						output_dv = output_dv + to_string(m[nDims-1][j]) + "," ;
				output_dv = output_dv + ".dat";
				writeSobols(output_dv,sobols);

				output_discrepancy = "data_discrepancy/" + to_string(nDims) + "D/optim/" + to_string(kmax) + "_" + to_string(tstDim) + "_";
				for(int j = 0; j < s[nDims-1]; j++)
					if(j == s[nDims-1] - 1)
						output_discrepancy = output_discrepancy + to_string(m[nDims-1][j]) ;
					else
						output_discrepancy = output_discrepancy + to_string(m[nDims-1][j]) + "," ;
				output_discrepancy = output_discrepancy + ".dat";
				writeStats(output_discrepancy, stats);
				found++;
				if(dbg_flag) cout  << itrial << " -> " << kmax << " -> " << found << endl;

				cout << ltm->tm_hour<<":"<<ltm->tm_min<<":"<<ltm->tm_sec << " npts = " << npts << " >>>>>>> output_dv into " << output_dv << endl;
			}
		}
	}

 }
