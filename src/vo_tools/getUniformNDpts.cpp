// getUniformNDpts.cpp
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

#include "../Samplers/SobolCascade.h"
#include "../Samplers/BlueTile.h"

using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)

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

uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216,33554432,67108864,134217728,268435456,536870912,1073741824,2147483648};
double d_powers_of_two[32] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216,33554432,67108864,134217728,268435456,536870912,1073741824,2147483648};
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
	unsigned int nDims = 6, firstDim = 0, npts = 16;
	std::string output_fname = "out.dat";
	std::string input_dir_vectors = "data/vo_sobol_init_tab.dat";
	bool dbg_flag = false, radomize_flag = true;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    uint32_t owen_tree_depth = 32, add_owen_tree_depth = 6;
    uint32_t nbits;
 	uint8_t NbThreadsMax = omp_get_max_threads();
	bool xy_flag = false;
	uint number_realizations=1;
	uint32_t seed = time(NULL);
	bool exportUTK=false;

	CLI::App app { "getUniformNDpts" };
	app.add_option("-d,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims)  );
	app.add_option("--xy,--xyPlusND", xy_flag, "this flag indicates whether XY is generated, default: "+ std::to_string(xy_flag)  );
	app.add_option("-n,--npts", npts, "number of points to generate, default: "+ std::to_string(npts) );
	app.add_option("-s,--seed", seed, "Random number generator seed. default: "+ std::to_string(seed)  );
	app.add_option("--id,--idv", input_dir_vectors, "input dir_vectors (ascii file), default: " + input_dir_vectors );
	app.add_option("-p,--permut", owen_permut_flag, "permut or not? when owen_permut_flag=0 ==> Sobol, default: "+ std::to_string(owen_permut_flag)  );
	app.add_option("-a,--add_owen_tree_depth", add_owen_tree_depth, "additional tree_depth wrt std owens depth, default: "+ std::to_string(add_owen_tree_depth) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)  );
  // Realizations and export "à la" utk
	app.add_flag("--exportUTK", exportUTK, "export as a utk point set (default: false)");
	app.add_option("-m", number_realizations, "number of realizations of the sampler (default: 1)");
	app.add_option("-o,--output", output_fname, "output fname, default: " + output_fname );
	CLI11_PARSE(app, argc, argv)

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

	//========================== allocation of tile_seeds
	std::mt19937_64 gen( seed );
    BlueTile seeds(TILE_SIZE, NDIMENSIONS, NSPP);
	//========================== end of allocation of tile_seeds

    if(xy_flag) nDims += 2;

	std::ofstream fout(output_fname, std::ofstream::out);
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
    if(dbg_flag) std::cout << "loading file " << input_dir_vectors << " " << d.size() << " entries." << std::endl;
    sobols.resize(m.size());

    double **points;
    points = (double**)malloc(sizeof(double*) * npts );
	for (unsigned int i = 0; i < npts; i++)
		points[i] = (double *)malloc(sizeof(double) * nDims );

    double kmax = 0;
    vector<t_stats> stats;
    init_sobols(sobols, d, s, a, m, seed, false);
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
    nbits = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));

	std::ofstream ofs(output_fname , std::ios_base::out);
	std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);  // uniform distribution of uint32_ts between 0 and 2^32-1
  
  //For all realizations
  for(int realization=0; realization < number_realizations; ++realization )
  {
    //REWRITING SEEDS
	seeds.random_init(gen);         // init par generateurs c++
    
    for (int ipt = 0; ipt < npts; ipt++)
		getUniformND(points[ipt] , sobols, seeds(0, 0, 0), nDims, ipt, nbits, owen_tree_depth, owen_permut_flag);
    
    if (exportUTK==false)
    {
      for (int ipt = 0; ipt < npts; ipt++) {
        for (unsigned int idim = 0; idim < nDims; idim++) {
          if(xy_flag)
            cout << std::setprecision(20) << std::fixed << points[ipt][(idim+nDims-2) % (nDims)] << " ";	// first: dims {nDims, nDims+1}, then {0, ..., nDims-1}
          else
            cout << std::setprecision(20) << std::fixed << points[ipt][idim] << " ";	// dims {0, ..., nDims-1}
          
        }
        cout << endl;
      }
      if (realization+1 != number_realizations)
        cout<<"#"<<std::endl;
    }
    else
    {
      //export utk
      for (int ipt = 0; ipt < npts; ipt++) {
        for (unsigned int idim = 0; idim < nDims; idim++) {
          if(xy_flag)
            ofs << std::setprecision(20) << std::fixed << points[ipt][(idim+nDims-2) % (nDims)] << " ";  // first: dims {nDims, nDims+1}, then {0, ..., nDims-1}
          else
            ofs << std::setprecision(20) << std::fixed << points[ipt][idim] << " ";  // dims {0, ..., nDims-1}
          
        }
        ofs << endl;
      }
      if (realization+1 != number_realizations)
        ofs<<"#"<<std::endl;
    }
  }
  ofs.close();
  
}

