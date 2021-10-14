#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
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
#define INTEGRAND_TYPE_HEAVISIDE	1
#define INTEGRAND_TYPE_GAUSS		2
//--------------------------------- structures
typedef uint32_t sobolInt;

//--------------------------------- global variables
double MSElimit_Gauss = 4.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];
uint32_t powers_of_four[16];

bool checkStratificationOne5D(
        const std::vector<SobolGenerator1D<uint32_t> >& std_sobols,
        const std::vector<SobolGenerator1D<uint32_t> >& nonstd_sobols,
		const unsigned int dim1,
		const unsigned int dim2,
		const unsigned int dim3,
		const unsigned int dim4,
		const unsigned int dim5,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
		const bool dbg_flag
		) {
	uint32_t x, y, z, w, u;
	uint32_t shift;
	uint32_t npts1D;
	unsigned char volumes[16][16][16][16][16];			// [16][16][16][16][16] for max octave = 20

	for (uint iOctave = octaveMin; iOctave <= octaveMax; iOctave++) {	// almost useless
		uint npts = powers_of_two[iOctave];
		switch (iOctave % 5) {
		case 0 :
			npts1D = powers_of_two[(iOctave) / 5];
			shift = (N_SIGNIFICANT_BITS - (iOctave) / 5);
			// init volumes
			for (uint iy = 0; iy < npts1D; iy++)
				for (uint ix = 0; ix < npts1D; ix++)
					for (uint iz = 0; iz < npts1D; iz++)
						for (uint iw = 0; iw < npts1D; iw++)
							for (uint iu = 0; iu < npts1D; iu++)
								volumes[iu][iw][iz][iy][ix] = 0;
			for (uint ipt = 0; ipt < npts; ipt++) {
				x = nonstd_sobols[dim1].getSobolInt(ipt) >> shift;
				y = nonstd_sobols[dim2].getSobolInt(ipt) >> shift;
				z = nonstd_sobols[dim3].getSobolInt(ipt) >> shift;
				w = nonstd_sobols[dim4].getSobolInt(ipt) >> shift;
				u = std_sobols[dim5].getSobolInt(ipt) >> shift;
				if (volumes[u][w][z][y][x] == 1) {
					return false;
				}
				volumes[u][w][z][y][x] = 1;
			}
			break;
		}
	}
	return true;
}	// checkStratificationOne5D

bool checkIntegrationGauss5D(
        const std::vector<SobolGenerator1D<sobolInt> >& std_sobols,
        const std::vector<SobolGenerator1D<sobolInt> >& nonstd_sobols,
        const unsigned int dim1,
        const unsigned int dim2,
        const unsigned int dim3,
        const unsigned int dim4,
        const unsigned int dim5,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        const unsigned int octaveStep,
        const unsigned int nTrials_integration,
		RNG& gen,
		const bool owen_permut_flag,
		const bool max_tree_depth_32_flag,
        const bool dbg_flag
		) {
	double ref_mse_Gauss[21] = {0.030328099227237725888, 0.0081343150413064422427, 0.0024911294750952087465, 0.0010547765401616429528, 0.00024353028376026040038, 4.9492811385256664514e-05, 1.6219815279860151115e-05,
			2.0662562335908198487e-06, 4.8885396746641327e-07, 7.66980820284343097e-08, 1.7505516524935886761e-08, 4.4308125348659952381e-09, 6.6400801035377593718e-10, 2.121309907194543851e-10,
			1.8375397748419005566e-11, 3.5980822865265122812e-12, 6.3422509626807047452e-13, 9.8925866002624597669e-14, 1.4054282739744241404e-14, 2.8270542897668849693e-15, 3.6932045931294861542e-16};
    vector<double> integral_estimations_Gauss;
    vector<t_pointND> mse_data_Gauss;
    double mse_accumulator_Gauss;
    unsigned int nTrials_integration_significant = nTrials_integration;;
	double kmax_Gauss = 0.;
	if (dbg_flag)
#pragma omp critical (output)
		cout << "Entering integration test for octaves from " << octaveMin << " to " << octaveMax << endl;
	for (int ioctave = octaveMin; ioctave <= octaveMax; ioctave += octaveStep) {
		size_t npts = powers_of_two[ioctave];
		uint32_t owen_tree_depth = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
		if(max_tree_depth_32_flag) owen_tree_depth = OWENPLUS_TREE_DEPTH;

		// init integral_estimations_Gauss
		integral_estimations_Gauss.clear();
        vector<VecX<5>> points(npts);
		for (int iset = 0; iset < nTrials_integration; iset++) {
            uint32_t seed1 = gen();
            uint32_t seed2 = gen();
            uint32_t seed3 = gen();
            uint32_t seed4 = gen();
            uint32_t seed5 = gen();
    		// prep pints
			for (sobolInt ipt = 0; ipt < npts; ipt++) {
				sobolInt IDcode;
				// dim dim1
				IDcode = nonstd_sobols[dim1].getSobolInt(ipt);
				if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed1, owen_tree_depth);
				points[ipt][0] = ((double) IDcode / (double) UINT32SOBOLNORM);
				// dim dim2
                IDcode = nonstd_sobols[dim2].getSobolInt(ipt);
                if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed2, owen_tree_depth);
                points[ipt][1] = ((double) IDcode / (double) UINT32SOBOLNORM);

				// dim dim3
				IDcode = nonstd_sobols[dim3].getSobolInt(ipt);
				if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed3, owen_tree_depth);
				points[ipt][2] = ((double) IDcode / (double) UINT32SOBOLNORM);
				// dim dim4
                IDcode = nonstd_sobols[dim4].getSobolInt(ipt);
                if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed4, owen_tree_depth);
                points[ipt][3] = ((double) IDcode / (double) UINT32SOBOLNORM);

				// dim dim5
                IDcode = std_sobols[dim5].getSobolInt(ipt);
                if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed5, owen_tree_depth);
                points[ipt][4] = ((double) IDcode / (double) UINT32SOBOLNORM);

			}
			// calculate integral_estimations
			double this_mse_Gauss = calculate_mse(points, INTEGRAND_TYPE_GAUSS);		// in SobolPlusPlus/src/Integration/Integration.h
	        integral_estimations_Gauss.push_back(this_mse_Gauss);

		}
		mse_accumulator_Gauss = 0.;
	    for (int i = 0; i < nTrials_integration_significant ; i++) {
	    	mse_accumulator_Gauss += integral_estimations_Gauss[i];
	    }
	    if(dbg_flag) {
	    	for (int i = 0; i < nTrials_integration_significant ; i++) cout << " " << integral_estimations_Gauss[i] << " " ;
	    	cout << endl;
	    }
		double mean_mse_Gauss = mse_accumulator_Gauss/(double)nTrials_integration_significant;
		double k_Gauss = mean_mse_Gauss / ref_mse_Gauss[ioctave];
		kmax_Gauss = max(kmax_Gauss, k_Gauss);

		if (dbg_flag) {
			cout  << dim1 << "," << dim2  << "," << dim3  << "," << dim4  << "," << dim5  << " =================>>> " << points.size() << " pts <<<================= "
				 << nTrials_integration << " trials -> " << " \t k_Gauss=" << k_Gauss << " \t kmax_Gauss=" << kmax_Gauss << endl;
		}
		if (kmax_Gauss > MSElimit_Gauss) {
			return false;
		}
	}

	if ( kmax_Gauss < MSElimit_Gauss) {
		return true;
	} else {
		return false;
	}
} // checkIntegrationGauss5D

bool checkStratificationAndIntegration5D(
		uint itrial,
        std::vector<SobolGenerator1D<sobolInt> >& std_sobols,
        std::vector<SobolGenerator1D<sobolInt> >& nonstd_sobols,
		const unsigned int dim1,
        const unsigned int dim2,
		const unsigned int dim3,
        const unsigned int dim4,
        const unsigned int dim5,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        std::ofstream& out,
        std::ofstream& out_strat,
		mt19937_64& gen,
        const unsigned int nTrials_integration,
		const bool owen_permut_flag,
		const bool max_tree_depth_32_flag,
        const bool dbg_flag) {

    bool test_strat5D = checkStratificationOne5D(std_sobols, nonstd_sobols, dim1, dim2, dim3, dim4, dim5, octaveMin, octaveMax, dbg_flag);
//        if(dbg_flag) cout << "===> " << itrial << " : " << dim1 << " " << dim2 << " " << dim3 << " " << dim4  << " -> " << test_strat5D << endl;
    if (test_strat5D) {
		cout << "===> "  << itrial << " : " << dim1 << " " << dim2 << " " << dim3 << " " << dim4 << " " << dim5
		     << "===> "  << itrial << " : " << nonstd_sobols[dim1].d << " " << nonstd_sobols[dim2].d << " " << nonstd_sobols[dim3].d << " " << nonstd_sobols[dim4].d << " " << std_sobols[dim5].d  << endl;

        std::uniform_int_distribution<uint32_t> unifFull(0);
        RNG randomGen;
        randomGen.seed(unifFull(gen));

        bool test_integration = false;
		test_integration = checkIntegrationGauss5D(std_sobols, nonstd_sobols, dim1, dim2, dim3, dim4, dim5, octaveMin, octaveMax, 1, nTrials_integration, randomGen,
    			owen_permut_flag, max_tree_depth_32_flag, dbg_flag);
    	if (test_integration) {
			out << nonstd_sobols[dim1] << endl;
			out << nonstd_sobols[dim2] << endl;
			out << nonstd_sobols[dim3] << endl;
			out << nonstd_sobols[dim4] << endl;
			out << std_sobols[dim5] << endl;
			out.flush();
			cout << " ############################################################################ " << endl;
			cout << nonstd_sobols[dim1] << endl;
			cout << nonstd_sobols[dim2] << endl;
			cout << nonstd_sobols[dim3] << endl;
			cout << nonstd_sobols[dim4] << endl;
			cout << std_sobols[dim5] << endl;
			cout << endl << " ############################################################################ " << endl;
			cout.flush();
			return true;
		}
    }
    return false;
}	//checkStratificationAndIntegration5D

int main(int argc, char **argv) {
	srand48( time(NULL) );
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i) powers_of_two[i] = 1u << i;
	int dim1, dim2, dim3, dim4, dim5;
	unsigned int nTrials_integration = 32;
	std::string filename = "out.dat", filename_strat = "out_strat.dat";
    uint8_t NbThreadsMax = omp_get_max_threads();
	bool dbg_flag = false;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    bool all_test_are_OK;
	std::string dir_vectors_fname = "data/data_SobolPlusPlus/2D+2D/best_found_2D+2D_upto12.dat";
	std::string std_dir_vectors_fname = "data/sobol_init_tab.dat";
	unsigned int octaveMin = 2, octaveMax = 8;	// up to 16K spp only

//	MSElimit_Gauss = 2.25;	// upto octaveMax = 10
	MSElimit_Gauss = 1.5;	// upto octaveMax = 8

	uint32_t seed = time(NULL); // 13374269;
	CLI::App app { "checkerStrat2D" };
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename );
//	app.add_option("--out_strat", filename_strat, "output filename_strat (ascii file), default: " + filename_strat );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: " + std::to_string(seed) );
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );
	app.add_option("--std_dir_vectors", std_dir_vectors_fname, "std_dir_vectors_fname (ascii file), default: " + std_dir_vectors_fname );
	app.add_option("--octaveMin", octaveMin, "min octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMin) );
	app.add_option("--octaveMax", octaveMax, "max octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMax) );
	app.add_option("--nTrials_integration", nTrials_integration, "nTrials_integration, default: " + std::to_string(nTrials_integration) );
	app.add_option("-p,--permut", owen_permut_flag, "owen_permut_flag, default: " + std::to_string(owen_permut_flag) );
	app.add_option("-m,--max_tree_depth_32", max_tree_depth_32_flag, "Owen or OwenPlus ? default: OwenPlus max_tree_depth_32_flag="+ std::to_string(max_tree_depth_32_flag)+"" );
	app.add_option("--MSElimit_Gauss", MSElimit_Gauss, "MSE goodness threshold during integration, default: " + std::to_string(MSElimit_Gauss) );
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + std::to_string(NbThreadsMax) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: " + std::to_string(dbg_flag) );
	CLI11_PARSE(app, argc, argv)

    std::uniform_int_distribution<uint32_t> unif01(0, 1);

	srand48(seed);
	std::ofstream out(filename, std::ofstream::app);
	std::ofstream out_strat;

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

    cout << "effective number of threads = "<< int(NbThreadsMax) << endl;
    cout << "octaveMin = "<< int(octaveMin) << "  " << "octaveMax = "<< int(octaveMax) << endl;
    cout << "MSElimit_Gauss = "<< MSElimit_Gauss << endl;

    std::vector<SobolGenerator1D<uint32_t> > std_sobols, nonstd_sobols;	// array of sobol data per dim

    cout << "Loading std Sobol dir_vectors from " << std_dir_vectors_fname << endl;
    loadSobolsFromFile(std_dir_vectors_fname, std_sobols,seed);		// read std_sobols from file
    cout << "Loading std Sobol dir_vectors from " << std_dir_vectors_fname << " of size " << std_sobols.size() << endl;

    cout << "Loading non-std Sobol dir_vectors from " << dir_vectors_fname << endl;
    loadSobolsFromFile(dir_vectors_fname, nonstd_sobols,seed);		// read std_sobols from file
    cout << "Loading non-std Sobol dir_vectors from " << dir_vectors_fname << " of size " << nonstd_sobols.size() << endl;

    std::mt19937_64 gen( time(NULL) );

    int n4tuples = nonstd_sobols.size() / 4;
    uint itrial = 0;
    while (true) {
    	itrial++;
        dim1 =  4 * floor(n4tuples*drand48());
        dim2 = dim1 + 1;
        dim3 = dim1 + 2;
        dim4 = dim1 + 3;

        dim5 =  floor( std_sobols.size() * drand48() );
		unsigned int second_entry = unif01(gen);
		std_sobols[dim5].m[1] = (second_entry << 1) ^ 1;
		for (int k = 2; k < std_sobols[dim5].s; ++k) {
			//Random in [0, 2^k[
			std::uniform_int_distribution<uint32_t> unifFull(0, powers_of_two[k] - 1);
			unsigned int i = unifFull(gen);
			std_sobols[dim5].m[k] = (i << 1) ^ 1;
		}
		std_sobols[dim5].generateMatrix();
        checkStratificationAndIntegration5D( itrial, std_sobols, nonstd_sobols,
       		 dim1,dim2,dim3,dim4,dim5,octaveMin,octaveMax,
			 out,out_strat, gen,
			 nTrials_integration,owen_permut_flag,max_tree_depth_32_flag,dbg_flag);
	}
}
