// demoCascadedSobol.cpp
// vo April 2021

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <cassert>
#include "CLI11.hpp"

#include "../Integration/Integration.h"
//#include "../Points/VecX.h"

#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

#include "../Samplers/SobolCascade.h"
#include "../Samplers/BlueTile.h"

using namespace std;

//--------------------------------- constants
#define N_SIGNIFICANT_BITS CHAR_BIT*sizeof(sobolInt)

#define TILE_SIZE	32
#define NDIMENSIONS 4
#define NSPP 2

//--------------------------------- structures
typedef uint32_t sobolInt;

struct type_pointND {double * coord ; };
struct type_pointND_int {uint32_t * coord ; };

//--------------------------------- global variables
struct t_sobol {
    int d, s, a;
    std::array<uint32_t, CHAR_BIT*sizeof(uint32_t)> m;
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
	    const bool dbg_flag = false ) 
{
    std::ofstream out_file(out_discrepancy_fname);
    if (!out_file.is_open()) {
        std::cerr << out_discrepancy_fname << " cannot open." << std::endl;
        exit(1);
    }
	out_file << "#Nbpts	#Mean	#Var	#Min	#Max	#Analytical	#MSE	#NbPtsets	#Nbintegrands\n" ;
	for (int i = 0; i < stats.size(); i++) {
		out_file << stats[i].Nbpts << " " << stats[i].Mean << " " << stats[i].Var << " " << stats[i].Min << " "
			<< stats[i].Max << " " << stats[i].Analytical << " " << stats[i].MSE << " " << stats[i].NbPtsets << " " << stats[i].k << endl;
	}
    if (dbg_flag) std::cout << "writing into file " << out_discrepancy_fname << " " << stats.size() << " entries." << std::endl;
	out_file.close();
}

inline 
void getOwenND_with_seed( VecX<NDIMENSIONS>& pt_ND,
    const std::vector<SobolGenerator1D<uint32_t> >& sobols,
    const uint32_t *seeds,
    const uint32_t nDims,
    const uint32_t n,
    const uint32_t nbits,
    const uint32_t owen_tree_depth = 32,
    const bool owen_permut_flag = true,
    const uint32_t first_dim = 2) 
{
    assert(first_dim + nDims <= sobols.size());
	for (unsigned int idim = 0; idim < nDims; idim++) {
	    uint32_t IDcode = n;		// dim 0: take n-th point as index -> into 32-bit integer IDcode
		IDcode = sobols[first_dim + idim].getSobolInt(IDcode);	// radix-inversion + sobol
		uint32_t res_IDcode = IDcode;				// used to calculate the resulting value
	    if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
        	res_IDcode = OwenScrambling(res_IDcode, seeds[idim], owen_tree_depth);
	    pt_ND[idim] = ((double) res_IDcode / (double) UINT32SOBOLNORM);	// final value (double) for this dimension
	}
}	// getUniformND


int main(int argc, char **argv) 
{
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i) 
        powers_of_two[i] = 1u << i;
    
	unsigned int nDims = NDIMENSIONS, firstDim = 0, npts, npts_max = 256;
	
    std::string output_fname = "out.dat";
	std::string input_cascaded_dir_vectors = "data/vo_sobol_init_tab.dat";
	std::string input_dir_vectors = "data/sobol_init_tab.dat";
    
	bool dbg_flag = false, radomize_flag = true;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    uint32_t owen_tree_depth = 32, add_owen_tree_depth = 6;
    uint32_t nbits;
 	uint8_t NbThreadsMax = omp_get_max_threads();
	bool xy_flag = false;		// when true, additional coordinates x y are generated before asked nDims
	uint32_t seed = time(NULL);
    vector<t_stats> stats;
    int integrandType = 1; // gaussians

	CLI::App app { "demoCascadedSobol" };
//	app.add_option("-d,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims)  );
	app.add_option("--xy,--xyPlusND", xy_flag, "this flag indicates whether XY is generated, default: "+ std::to_string(xy_flag)  );
	app.add_option("-n,--npts_max", npts_max, "number of points to generate, default: "+ std::to_string(npts_max) );
	app.add_option("-s,--seed", seed, "Random number generator seed. default: "+ std::to_string(seed)  );
//	app.add_option("-o,--output", output_fname, "output fname, default: " + output_fname );
	app.add_option("--cdirs", input_cascaded_dir_vectors, "cascaded input dir_vectors (ascii file), default: " + input_cascaded_dir_vectors );
	app.add_option("--dirs", input_dir_vectors, "input dir_vectors (ascii file), default: " + input_dir_vectors );
    
	app.add_option("-p,--permut", owen_permut_flag, "permut or not? when owen_permut_flag=0 ==> Sobol, default: "+ std::to_string(owen_permut_flag)  );
//	app.add_option("-a,--add_owen_tree_depth", add_owen_tree_depth, "additional tree_depth wrt std owens depth, default: "+ std::to_string(add_owen_tree_depth) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)  );
	CLI11_PARSE(app, argc, argv)

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

	if(xy_flag) nDims += 2;


	//============================= CascadedSobolTab
	//========================== initialization of tile_seeds
	std::mt19937_64 gen( seed );
    
    BlueTile seeds(TILE_SIZE, NDIMENSIONS, NSPP);
    //~ seeds.random_init(seed);    // init par RNG de Random.h
    seeds.random_init(gen);         // init par generateurs c++
    seeds.write("random_seeds.dat");
    
    
	//========================== loading direction vectors
	std::ofstream fout(output_fname, std::ofstream::out);
	std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    {
        sobols.clear();
        loadSobolsFromFile(input_cascaded_dir_vectors, sobols, dbg_flag);		// read cascaded sobols from file and fill appropriate structures
        
        if(dbg_flag)
            for(int i = 0; i < nDims; i++) {
                cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
                for(int j = 0; j < sobols[i].s; j++) cout << " " << sobols[i].m[j];
                cout  << endl;
            }
    }
    
	//========================== making pointsets and stats
    vector<VecX<NDIMENSIONS>> points;//(npts_max);
    
    FILE *mse_output_cascade= fopen("csobol.dat", "wt");
    assert(mse_output_cascade);
    
	for (unsigned int nbits = 2; nbits <= 12; nbits += 2) {
		npts = powers_of_two[nbits];
//		if(owen_permut_flag) owen_tree_depth = nbits + add_owen_tree_depth;		// owen_tree_depth is fixed to 32
        
		double mse_cumul = 0.;
		for(int iy = 0; iy < seeds.height(); iy++) 
        {
			for(int ix = 0; ix < seeds.width(); ix++) 
            {
				points.resize(npts);
				for (int ipt = 0; ipt < npts ; ipt++) {
                    //~ seeds(ix, iy, id, is) = 12345;      // modifie le seed pour le pixel ix, iy, dimensions id, sample is
                    //~ uint32_t seed= seeds(ix, iy, id, is);   // recupere le seed
                    //~ const uint32_t *seed= seeds(ix, iy, 0); // renvoie le tableau de seed 
                    
					getUniformND(points[ipt].data(), sobols, seeds(ix, iy, 0), nDims, ipt, nbits, owen_tree_depth, owen_permut_flag);
				}
				double mse = calculate_mse(points, integrandType, 1024, 13374269); // deterministic selection
				mse_cumul += mse;
                
				if (dbg_flag) {
					for (int ipt = 0; ipt < points.size() ; ipt++) {
						for (int idim = 0; idim < NDIMENSIONS ; idim++)
							cout << points[ipt][idim] << " " ;
						cout << endl;
					}
					cout << nDims << "D -> " << npts << " " << ix << " " << iy << " " << " -> " << mse << endl;
				}
			}
		}
        mse_cumul /= (double)(TILE_SIZE*TILE_SIZE);
        
		cout << npts << " " << mse_cumul << endl;
        fprintf(mse_output_cascade, "%d %.10f\n", npts, mse_cumul);
	}
    
    fclose(mse_output_cascade);
	cout<< endl;
    
    
	//============================= the same, but with owen
	//========================== loading direction vectors
    {
        sobols.clear();    
        loadSobolsFromFile(input_dir_vectors, sobols, dbg_flag);		// read sobols from file and fill appropriate structures
        
        if(dbg_flag)
            for(int i = 0; i < nDims; i++) {
                cout << i << " : " << sobols[i].d << " " << sobols[i].s << " " << sobols[i].a << " ";
                for(int j = 0; j < sobols[i].s; j++) cout << " " << sobols[i].m[j];
                cout  << endl;
            }
    }
    
	//========================== making pointsets and stats

    FILE *mse_output_owen= fopen("owen.dat", "wt");
    assert(mse_output_owen);
    
	for (unsigned int nbits = 2; nbits <= 12; nbits += 2) {
		npts = powers_of_two[nbits];
//		if(owen_permut_flag) owen_tree_depth = nbits;	// std sobol/owen
		owen_tree_depth = 32;
        
		double mse_cumul = 0.;
		for(int iy = 0; iy < seeds.height(); iy++) 
        {
			for(int ix = 0; ix < seeds.width(); ix++) 
            {
				points.resize(npts);
				for (int ipt = 0; ipt < npts ; ipt++) 
                    getOwenND_with_seed(points[ipt], sobols, seeds(ix, iy, 0), nDims, ipt, nbits, owen_tree_depth, owen_permut_flag);
                
				//~ double mse = calculate_mse(points, integrandType, 1024); // random selection
				double mse = calculate_mse(points, integrandType, 1024, 13374269); // deterministic selection
				mse_cumul += mse;
                
				if (dbg_flag) {
					for (int ipt = 0; ipt < points.size() ; ipt++) {
						for (int idim = 0; idim < NDIMENSIONS ; idim++)
							cout << points[ipt][idim] << " " ;
						cout << endl;
					}
					cout << nDims << "D -> " << npts << " " << ix << " " << iy << " " << " -> " << mse << endl;
				}
			}
		}
        mse_cumul /= (double)(TILE_SIZE*TILE_SIZE);
        
        cout << npts << " " << mse_cumul << endl;
        fprintf(mse_output_owen, "%d %.10f\n", npts, mse_cumul);
	}
    
    fclose(mse_output_owen);
	cout<< endl;

}
