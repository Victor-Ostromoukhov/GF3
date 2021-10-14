// getUniformNDdiscrepancy.cpp
// vo March 2021

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include "CLI11.hpp"

#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)

//--------------------------------- structures
typedef uint32_t sobolInt;

struct type_pointND {double * coord ; };

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
inline double calculate_GeneralizedL2Discrepancy(const unsigned int nDims, const std::vector<type_pointND>& points )
{
	uint N = points.size();
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
			double uij = points[i].coord[j];
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
			double uij = points[i].coord[j];
			double uiprimej = points[iprime].coord[j];
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

int main(int argc, char **argv) {
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i) powers_of_two[i] = 1u << i;
	CLI::App app { "owen" };
	unsigned int nDims = 6, firstDim = 0, npts = 16;
	std::string output_discrepancy;
	std::string input_dir_vectors = "data/vo_sobol_init_tab.dat";
	bool dbg_flag = true;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    uint32_t set_owen_tree_depth = 0;
    uint32_t nbits;
    uint32_t ocatve_from = 4, octave_to = 20, octave_step = 2, dv_length = -1;
	uint8_t NbThreadsMax = omp_get_max_threads();

	uint32_t seed = time(NULL);
	app.add_option("-d,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims)+")" );
	app.add_option("--id,--idv", input_dir_vectors, "input dir_vectors (ascii file), default: " + input_dir_vectors );
    app.add_option("--of,--ocatve_from", ocatve_from, "starting octave, default: "+ std::to_string(ocatve_from)+")" );
    app.add_option("--ot,--ocatve_to", octave_to, "ending octave, default: "+ std::to_string(octave_to)+")" );
    app.add_option("--os,--octave_step", octave_step, "octave step, default: "+ std::to_string(octave_step)+")" );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)+")" );
	CLI11_PARSE(app, argc, argv)

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

    srand48(seed);
    std::mt19937_64 gen(seed);
	vector<type_pointND> points;

//	std::ofstream fout(input_dir_vectors, std::ofstream::out);
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
    std::cout << "loading file " << input_dir_vectors << " " << d.size() << " entries." << std::endl;
    sobols.resize(m.size());

    double kmax = 0;
    vector<t_stats> stats;
    init_sobols(sobols, d, s, a, m, seed, false);
    if(dbg_flag)
	    for(int i = 0; i < nDims; i++) {
	        cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
	        for(int j = 0; j < sobols[i].d; j++) cout << " " << sobols[i].m[j];
	        cout  << endl;
	    }
	for (int octave = ocatve_from; octave <= octave_to; octave += octave_step) {
	    struct t_stats this_stats;
		npts = powers_of_two[octave];
		nbits = octave;
		points.clear();
		for (int ipt = 0; ipt < npts; ipt++) {
			type_pointND pt_ND;
			pt_ND.coord = (double *)malloc(sizeof(double) * nDims);
			uint ind = ipt;
			for (unsigned int idim = 0; idim < nDims; idim++) {
				ind = (sobols[idim].getSobolInt(ind) >> (32-nbits));
				pt_ND.coord[idim] = ((double) ind  / (double) npts);
			}
			points.push_back( pt_ND );
		}
		double discrepancy = calculate_GeneralizedL2Discrepancy(nDims, points);
		double k;
		k = discrepancy / discrepancy_refs[nDims][octave];
		kmax = max(k, kmax);
		if(dbg_flag)
			cout << "npts = " << npts << " \t d -> " << discrepancy << " \t kmax -> " << kmax  << " \tk -> " << k << endl;
		this_stats.Nbpts = npts;
		this_stats.Mean = discrepancy;
		this_stats.Var	=	this_stats.Min	=	this_stats.Max	=	this_stats.Analytical	=	this_stats.MSE	=	0;
		this_stats.NbPtsets	= 1;
		this_stats.k = k;
		stats.push_back ( this_stats );
	}
	output_discrepancy = "data_discrepancy/" + to_string(nDims) + "D/UniformND_" + to_string(kmax) + ".dat";
	cout << ">>>>>>> output_discrepancy into " << output_discrepancy << endl;
	writeStats(output_discrepancy, stats);
 }

//cout << std::bitset<32>(IDcode ) << " \t";

