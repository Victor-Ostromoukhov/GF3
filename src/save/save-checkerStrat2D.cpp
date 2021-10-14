#include <iostream>
#include <string>
#include <fstream>
#include <omp.h>


#include "CLI11.hpp"
#include "sobol+owen.hpp"

using namespace std;

//--------------------------------- constants
#define NDIMS 2

//--------------------------------- structures
struct t_pointND {double coord[NDIMS] ; };

struct t_GaussianStruct2D {double integral; double mu[2] ; double mxCInv[2 * 2] ; };
struct t_GaussianStruct3D {double integral; double mu[3] ; double mxCInv[3 * 3] ; };
struct t_GaussianStruct4D {double integral; double mu[4] ; double mxCInv[4 * 4] ; };
struct t_GaussianStruct6D {double integral; double mu[6] ; double mxCInv[6 * 6] ; };

//--------------------------------- global variables
//uint32_t powers_of_two[32]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216,33554432,67108864,134217728,268435456,536870912,1073741824,2147483648};
// powers_of_two are aleady declared in "sobol+owen.hpp"
unsigned char squares[1024][1024];		// for static allocation 2^20 max
uint32_t owen_tree_depth = 32;
double MSElimit = 3.;
#include "includes/GaussianStruct2D.hpp"
#include "includes/GaussianStruct3D.hpp"
#include "includes/GaussianStruct4D.hpp"
#include "includes/GaussianStruct6D.hpp"

//--------------------------------- routines related to integration
double getMultivariateGaussian(const int nDims, const double *pt_ND, const double *mu, double *mxCInv) {
    double accumulator = 0.;
	for (int row = 0; row < nDims; row++) {
		for (int col = 0; col < nDims; col++)accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
	}
	return exp(-.5 * accumulator);
}	// getMultivariateGaussian


double calculate_mse_Gaussian(const int nDims, vector<t_pointND>& points, int npts) {
    double mse_accumulator = 0.;
    int nintegrands=1024;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(16); // Use 16 threads for all consecutive parallel regions#pragma omp parallel for schedule(dynamic) reduction(+:mse_accumulator)
#pragma omp parallel reduction(+:mse_accumulator)
//#pragma omp parallel for schedule(dynamic) reduction(+:mse_accumulator)
    for (int iintegrands = 0; iintegrands < nintegrands; iintegrands ++) {
        double integral, *mu, *mxCInv;
        switch (nDims) {
        case 2 : integral = tab_GaussianStruct2D[iintegrands].integral;
    		mu = tab_GaussianStruct2D[iintegrands].mu;
    		mxCInv = tab_GaussianStruct2D[iintegrands].mxCInv;
        	break;
        case 3 : integral = tab_GaussianStruct3D[iintegrands].integral;
        	mu = tab_GaussianStruct3D[iintegrands].mu;
    		mxCInv = tab_GaussianStruct3D[iintegrands].mxCInv;
        	break;
        case 4 : integral = tab_GaussianStruct4D[iintegrands].integral;
			mu = tab_GaussianStruct4D[iintegrands].mu;
			mxCInv = tab_GaussianStruct4D[iintegrands].mxCInv;
			break;
        case 6 : integral = tab_GaussianStruct6D[iintegrands].integral;
			mu = tab_GaussianStruct6D[iintegrands].mu;
			mxCInv = tab_GaussianStruct6D[iintegrands].mxCInv;
			break;
        default: cout << "not implemented" << endl;
        	exit(1);
            break;
        }
        double accumulator = 0.;
        for (int ipt = 0; ipt < points.size() ; ipt++) {
			accumulator += getMultivariateGaussian(nDims, points[ipt].coord, mu, mxCInv);
        }
        accumulator /= (double)npts;
        double mse = (integral-accumulator)*(integral-accumulator);
//        integral_estimations.push_back(accumulator);
        mse_accumulator += mse;
    }
    return mse_accumulator / (double)(nintegrands);
} // calculate_mse_Gaussian

bool checkIntegrationOnePairOfDims(
		t_sampler_env *sampler_env,
		const unsigned int dim1,
		const unsigned int dim2,
		const unsigned int octaveMax,
		const unsigned int octaveStep,
		const bool dbg_flag
		) {
	struct timeval tp; gettimeofday(&tp, NULL);
	uint global_random_seed = tp.tv_usec + getpid();
	srand48(global_random_seed);
	double mse_accumulator = 0.;
	int ntrials = 32;

	double ref_mse[15] = {
			0.082541622832531044706705,		// ocatve 0
			0.019356097022389708722789,		// ocatve 1
			0.003925665357197775529211,		// ocatve 2
			0.000720739293620302894205,		// ocatve 3
			0.0001215184018476366665997,	// ocatve 4
			0.00001921158804074820053581,	// ocatve 5
			2.92779851792596907458e-6,		// ocatve 6
			4.3154559275337646244638e-7,	// ocatve 7
			6.216779003101871597345e-8,		// ocatve 8
			8.81220066619763449804e-9,		// ocatve 9
			1.22708716714593846023e-9,		// ocatve 10
			1.696014440334198589017e-10,	// ocatve 11
			2.31838467946709772845e-11,		// ocatve 12
			3.1329021425225278940123e-12,	// ocatve 13
			4.2443156384433291527584e-13};	// ocatve 14

	vector<t_pointND> points;
	t_pointND pt_ND;

	double kmax = 0.;
	for (int ioctave = 2; ioctave <= octaveMax; ioctave += octaveStep) {
		int npts = powers_of_two[ioctave];
		mse_accumulator = 0.;
		for (int iset = 0; iset < ntrials; iset++) {
			points.clear();
			sampler_init_seeds_one_entry(sampler_env, 7531*(iset+1), dim1, dbg_flag);
			sampler_init_seeds_one_entry(sampler_env, 5731*(iset+1), dim2, dbg_flag);
			for (int ipt = 0; ipt < npts; ipt++) {
				uint32_t IDcode;
				// dim dim1
				IDcode = get_1pt_sobol1D(dim1, ipt);
				IDcode = getIDcodePermutedBinary(sampler_env, IDcode, dim1, owen_tree_depth);
				pt_ND.coord[0] = ((double) IDcode / (double) UINT32NORM);
				// dim dim2
				IDcode = get_1pt_sobol1D(dim2, ipt);
				IDcode = getIDcodePermutedBinary(sampler_env, IDcode, dim2, owen_tree_depth);
				pt_ND.coord[1] = ((double) IDcode / (double) UINT32NORM);
				points.push_back( pt_ND );
			}
			if (dbg_flag) {
		        for (int ipt = 0; ipt < points.size() ; ipt++) {
		        	pt_ND = points[ipt];
					cout << floor(npts*pt_ND.coord[0]) << " " << floor(npts*pt_ND.coord[1]) << endl;
		        }
		        cout  << endl;
			}
			double this_mse = calculate_mse_Gaussian(2, points, npts);
			mse_accumulator += this_mse;
		}
		double mean_mse = mse_accumulator/(double)ntrials;
		double k = mean_mse / ref_mse[ioctave];
		kmax = max(kmax,k);
		if(dbg_flag)
			cout << dim1 << "," << dim2 << " " << points.size() << " -> " << " " << mean_mse << " \tk=" << k  << " \tkmax=" << kmax << "/" << MSElimit << " -> " << (kmax > MSElimit) << endl;
		if (kmax > MSElimit) return false;
	}
	if (kmax < MSElimit) {
		return true;
	} else {
		return false;
	}
} // checkIntegrationOnePairOfDims

//--------------------------------- routines related to checkStratification
bool checkStratificationOnePairOfDims(
		const t_sampler_env *sampler_env,
		const unsigned int dim1,
		const unsigned int dim2,
		const unsigned int octaveMax,
		const bool squares_only_flag,
		const bool dbg_flag
		) {
	bool found;
	uint32_t x, y, shift, npts1D;
	for (uint iOctave = 1; iOctave <= octaveMax; iOctave++) {
		uint npts = powers_of_two[iOctave];
		found = true;
		if ((iOctave & 1) == 0) { //enven octave
			npts1D = powers_of_two[(iOctave) / 2];
			shift = (N_SIGNIFICANT_BITS - (iOctave) / 2);
			// init squares
			for (uint iy = 0; iy < npts1D; iy++)
				for (uint ix = 0; ix < npts1D; ix++)
					squares[iy][ix] = 0;
			for (uint ipt = 0; ipt < npts; ipt++) {
				x = get_1pt_sobol1D(dim1, ipt) >> shift;
				y = get_1pt_sobol1D(dim2, ipt) >> shift;
				if (squares[y][x] == 1) {
					found = false;
					return found;
				}
				squares[y][x] = 1;
			}
			if (!found)
				break;
		} else { // odd octave
			if (!squares_only_flag) {
				npts1D = powers_of_two[(iOctave + 1) / 2];
				shift = (N_SIGNIFICANT_BITS - (iOctave + 1) / 2);
				// init squares
				for (uint iy = 0; iy < npts1D; iy++)
					for (uint ix = 0; ix < npts1D; ix++)
						squares[iy][ix] = 0;
				for (uint ipt = 0; ipt < npts; ipt++) {
					if (shift == N_SIGNIFICANT_BITS - 1) {
						x = get_1pt_sobol1D(dim1, ipt) >> shift;
						y = 0;
					} else {
						x = get_1pt_sobol1D(dim1, ipt) >> shift;
						y = get_1pt_sobol1D(dim2, ipt) >> (shift + 1);
					}
					if (squares[y][x] == 1) {
						found = false;
						return found;
					}
					squares[y][x] = 1;
				}
				// init once again
				for (uint iy = 0; iy < npts1D; iy++)
					for (uint ix = 0; ix < npts1D; ix++)
						squares[iy][ix] = 0;
				for (uint ipt = 0; ipt < npts; ipt++) {
					if (shift == N_SIGNIFICANT_BITS - 1) {
						x = 0;
						y = get_1pt_sobol1D(dim2, ipt) >> shift;
					} else {
						x = get_1pt_sobol1D(dim1, ipt) >> (shift + 1);
						y = get_1pt_sobol1D(dim2, ipt) >> shift;
					}
					if (squares[y][x] == 1) {
						found = false;
						return found;
					}
					squares[y][x] = 1;
				}
			}	// if (!squares_only_flag) {
			if (!found)
				break;
		}
	}
	return found;
}	// checkStratificationOnePairOfDims

void checkStratification(
		t_sampler_env *sampler_env,
		const unsigned int dim1,
		const unsigned int dim2,
		const unsigned int nTrialsdim1,
		const unsigned int nTrialsdim2,
		const unsigned int octaveMax,
		const std::string &filename,
		const bool dbg_flag, const bool squares_only_flag) {
	cout << "Checking Sobol indices (" << dim1 << "," << dim2 << ") for nTrialsdim1=" << nTrialsdim1 << " nTrialsdim2=" << nTrialsdim2 << " and npts between 2^1 and 2^" << octaveMax << " output into " << filename << endl;
	auto npairs = 0;
	std::ofstream outfile;
	outfile.open(filename, std::ios_base::app);
//#pragma parallel for schedule(dynamic)
	for (auto itrial1 = 1; itrial1 <= nTrialsdim1; itrial1++) {
		if(dbg_flag) cout << dim1 << "," << dim2 << "   \t" << itrial1 << "/" << nTrialsdim1 << endl;
		for (auto itrial2 = 1; itrial2 <= nTrialsdim2; itrial2++) {
			unsigned int second_entry = round( drand48() );
			sobol_mk[dim1][1] = (second_entry << 1) ^ 1;
			sobol_mk[dim2][1] = ((1-second_entry) << 1) ^ 1;
			for (int k = 2; k < sobol_sj[dim1]; ++k) {
				unsigned int i = floor((double)powers_of_two[k] * drand48());
				sobol_mk[dim1][k] = (i << 1) ^ 1;
			}
			for (int k = 2; k < sobol_sj[dim2]; ++k) {
				unsigned int i = floor((double)powers_of_two[k] * drand48());
				sobol_mk[dim2][k] = (i << 1) ^ 1;
			}
			if (checkStratificationOnePairOfDims(sampler_env, dim1, dim2, octaveMax, squares_only_flag, dbg_flag)) {
				//We have found a good seed par
//        #pragma omp critical
				cout << dim1 << "," << dim2 << "   \t" << itrial1 << "/" << nTrialsdim1 << "   \t" << itrial2 << "/" << nTrialsdim2 << endl;
				if(dbg_flag) {
					cout << sobol_dj[dim1] << "\t" << sobol_sj[dim1] << "\t" << sobol_aj[dim1]  << "\t";
					for (uint i = 0; i < sobol_sj[dim1]; ++i) cout << sobol_mk[dim1][i] << " ";
					cout << endl;
					cout << sobol_dj[dim2] << "\t" << sobol_sj[dim2] << "\t" << sobol_aj[dim2]  << "\t";
					for (uint i = 0; i < sobol_sj[dim2]; ++i) cout << sobol_mk[dim2][i] << " ";
					cout << endl;
				}
				if( checkIntegrationOnePairOfDims(sampler_env,dim1,dim2,octaveMax,1,dbg_flag) ) {	// check with octaveStep=1
					npairs++;
					cout << npairs << " ############################################################################ " << filename << endl;
					cout << sobol_dj[dim1] << "\t" << sobol_sj[dim1] << "\t" << sobol_aj[dim1]  << "\t";
					for (uint i = 0; i < sobol_sj[dim1]; ++i) cout << sobol_mk[dim1][i] << " ";
					cout << endl;
					cout << sobol_dj[dim2] << "\t" << sobol_sj[dim2] << "\t" << sobol_aj[dim2]  << "\t";
					for (uint i = 0; i < sobol_sj[dim2]; ++i) cout << sobol_mk[dim2][i] << " ";
					cout << endl;
					cout << npairs << " ############################################################################ " << filename << endl;
					outfile << sobol_dj[dim1] << "\t" << sobol_sj[dim1] << "\t" << sobol_aj[dim1]  << "\t";
					for (uint i = 0; i < sobol_sj[dim1]; ++i) outfile << sobol_mk[dim1][i] << " ";
					outfile << endl;
					outfile << sobol_dj[dim2] << "\t" << sobol_sj[dim2] << "\t" << sobol_aj[dim2]  << "\t";
					for (uint i = 0; i < sobol_sj[dim2]; ++i) outfile << sobol_mk[dim2][i] << " ";
					outfile << endl;
				}
			}
		}
	}
	std::cout << "found " << npairs << " pairs of good seeds" << std::endl;
	outfile.close();
} // checkStratification

int main(int argc, char **argv) {
	srand48( time(NULL) );
	CLI::App app { "checkerStrat2D" };
	unsigned int dim1=1, nTrialsdim1 = 8192, nTrialsdim2 = 8192; // 4096 8192 32768;
	unsigned int octaveMax = 14;
	std::string filename = "out.dat", dir_vectors_fname = "data/sobol_init_tab.dat";
	bool dbg_flag = false, squares_only_flag = false;
	app.add_option("-i,--dim1", dim1, "Index of the first Sobol pair (==dim1) convetion: dim1=0 : van der Corput")->required();
	app.add_option("-o,--output", filename, "output filename (ascii file), default: out.dat");
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default:data/sobol_init_tab.dat");
	app.add_option("--octaveMax", octaveMax, "max octave value (will check all 2^i spp for i {2..octaveMax}), default: 14");
	app.add_option("--nTrialsdim1", nTrialsdim1, "nTrialsdim1, default: 8192");
	app.add_option("--nTrialsdim2", nTrialsdim2, "nTrialsdim2, default: 8192");
	app.add_option("--squares_only", squares_only_flag, "squares_only_flag, default: false");
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: false");
	CLI11_PARSE(app, argc, argv);
	t_sampler_env *sampler_env = sampler_init(dir_vectors_fname, 0, false);

	for (unsigned int dim2 = dim1 + 1; dim2 <= dim1 + 64; ++dim2) {
		checkStratification(sampler_env, dim1, dim2, nTrialsdim1, nTrialsdim2, octaveMax, filename, dbg_flag, squares_only_flag);
	}

}
