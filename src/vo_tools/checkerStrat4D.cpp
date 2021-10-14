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
double MSElimit_Gauss = 5., kmax_best_found = 1000.;

uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];
uint32_t powers_of_four[16];

bool checkStratificationOne4D(
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
		const unsigned int dim1,
		const unsigned int dim2,
		const unsigned int dim3,
		const unsigned int dim4,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
		const bool dbg_flag
		) {
	uint32_t x, y, z, w;
	uint32_t shift;
	uint32_t npts1D;
	unsigned char volumes[32][32][32][32];			// [32][32][32][32] for max octave = 20

	for (uint iOctave = octaveMin; iOctave <= octaveMax; iOctave++) {	// almost useless
		uint npts = powers_of_two[iOctave];
		switch (iOctave % 4) {
		case 0 :
			npts1D = powers_of_two[(iOctave) / 4];
			shift = (N_SIGNIFICANT_BITS - (iOctave) / 4);
			// init volumes
			for (uint iy = 0; iy < npts1D; iy++)
				for (uint ix = 0; ix < npts1D; ix++)
					for (uint iz = 0; iz < npts1D; iz++)
						for (uint iw = 0; iw < npts1D; iw++)
							volumes[iw][iz][iy][ix] = 0;
			for (uint ipt = 0; ipt < npts; ipt++) {
				x = sobols[dim1].getSobolInt(ipt) >> shift;
				y = sobols[dim2].getSobolInt(ipt) >> shift;
				z = sobols[dim3].getSobolInt(ipt) >> shift;
				w = sobols[dim4].getSobolInt(ipt) >> shift;
				if (volumes[w][z][y][x] == 1) {
					return false;
				}
				volumes[w][z][y][x] = 1;
			}
			break;
//		case 1 :
//			npts1D = powers_of_two[(iOctave+1) / 4];
//			shift = (N_SIGNIFICANT_BITS - (iOctave+1) / 4);
//			// init volumes
//			for (uint iy = 0; iy < npts1D; iy++)
//				for (uint ix = 0; ix < npts1D; ix++)
//					for (uint iz = 0; iz < npts1D; iz++)
//						for (uint iw = 0; iw < npts1D; iw++)
//							volumes[iw][iz][iy][ix] = 0;
//			for (uint ipt = 0; ipt < npts; ipt++) {
//				x = (uint64_t)sobols[dim1].getSobolInt(ipt) >> (shift);
//				y = (uint64_t)sobols[dim2].getSobolInt(ipt) >> (shift+1);
//				z = (uint64_t)sobols[dim3].getSobolInt(ipt) >> (shift+1);
//				w = (uint64_t)sobols[dim4].getSobolInt(ipt) >> (shift+1);
////		cout <<npts<<" "<< x << y << z << w <<endl;
//				if (volumes[w][z][y][x] == 1) {
//					return false;
//				}
//				volumes[w][z][y][x] = 1;
//			}
//			// init volumes
//			break;
		}
	}
	return true;
}	// checkStratificationOne4D

bool checkIntegrationGauss4D(
        const std::vector<SobolGenerator1D<sobolInt> >& sobols,
        const unsigned int dim1,
        const unsigned int dim2,
        const unsigned int dim3,
        const unsigned int dim4,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        const unsigned int octaveStep,
        const unsigned int nTrials_integration,
		RNG& gen,
		const bool owen_permut_flag,
		const bool max_tree_depth_32_flag,
        const bool dbg_flag
		) {
	double ref_mse_Gauss[21] = {0.04295883754254002627, 0.0106874, 0.003053374311827227317, 0.0010457009493131707,
			0.0001813242382154007623, 0.00003528088862823913151, 4.504705933795904713e-6, 1.36306694581662749e-6,
			4.192419323186763693e-7, 4.542366570142065922e-8, 5.538604568594615886e-9, 1.0157253897795721308e-9,
			1.674155571890043e-10, 3.159625206020291656e-11, 4.959618654627772417e-12, 7.57747104199093188e-13,
			1.0934760577013432e-13, 1.811075245914993438e-14, 2.523103314999800797e-15, 3.486776547940874255e-16, 5.417049942566438627e-17};
    vector<double> integral_estimations_Gauss;
    vector<t_pointND> mse_data_Gauss;
    double mse_accumulator_Gauss;
    unsigned int nTrials_integration_significant = nTrials_integration;;
	double kmax_Gauss = 0.;
	if (dbg_flag)
#pragma omp critical (output)
		cout << "Entering integration test" << endl;
	for (int ioctave = octaveMin; ioctave <= octaveMax; ioctave += octaveStep) {
		size_t npts = powers_of_two[ioctave];
		uint32_t owen_tree_depth = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
		if(max_tree_depth_32_flag) owen_tree_depth = OWENPLUS_TREE_DEPTH;

		// init integral_estimations_Gauss
		integral_estimations_Gauss.clear();
        vector<VecX<4>> points(npts);
		for (int iset = 0; iset < nTrials_integration; iset++) {
            uint32_t seed1 = gen();
            uint32_t seed2 = gen();
            uint32_t seed3 = gen();
            uint32_t seed4 = gen();
    		// prep pints
			for (sobolInt ipt = 0; ipt < npts; ipt++) {
				sobolInt IDcode;
				// dim dim1
				IDcode = sobols[dim1].getSobolInt(ipt);
				if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed1, owen_tree_depth);
				points[ipt][0] = ((double) IDcode / (double) UINT32SOBOLNORM);
				// dim dim2
                IDcode = sobols[dim2].getSobolInt(ipt);
                if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed2, owen_tree_depth);
                points[ipt][1] = ((double) IDcode / (double) UINT32SOBOLNORM);

				// dim dim3
				IDcode = sobols[dim3].getSobolInt(ipt);
				if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed3, owen_tree_depth);
				points[ipt][2] = ((double) IDcode / (double) UINT32SOBOLNORM);
				// dim dim4
                IDcode = sobols[dim4].getSobolInt(ipt);
                if(owen_permut_flag) IDcode = OwenScrambling(IDcode, seed4, owen_tree_depth);
                points[ipt][3] = ((double) IDcode / (double) UINT32SOBOLNORM);

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
			cout  << dim1 << "," << dim2 << " / " << sobols.size() << " | " << sobols[dim1].d - 1 << "," << sobols[dim2].d - 1 << " =================>>> " << points.size() << " pts <<<================= "
				 << nTrials_integration << " trials -> " << " \t k_Gauss=" << k_Gauss << " \t kmax_Gauss=" << kmax_Gauss << endl;
		}
		if (kmax_Gauss > MSElimit_Gauss) {
			kmax_best_found = min(kmax_best_found, kmax_Gauss);
			return false;
		}
//		if(doubleCheck_flag) {
//			t_pointND this_measure;
//			this_measure.coord[0] = npts;
//			this_measure.coord[1] = mean_mse_Gauss;
//			mse_data_Gauss.push_back( this_measure );
//		}
	}

	kmax_best_found = min(kmax_best_found, kmax_Gauss);

	if ( kmax_Gauss < MSElimit_Gauss) {
//		if(doubleCheck_flag) {
//			string fname;
//			fname = "selected_dat_dv/dv_" + to_string(sobols[dim1].d) + "_" + to_string(sobols[dim1].s) + "_" + to_string(sobols[dim1].a);
//			for (uint i = 0; i < sobols[dim1].s; ++i) fname += ( "_" + to_string(sobols[dim1].m[i]) );
//			fname += ("," + to_string(sobols[dim2].d) + "_" + to_string(sobols[dim2].s) + "_" + to_string(sobols[dim2].a) );
//			for (uint i = 0; i < sobols[dim2].s; ++i) fname += ( "_" + to_string(sobols[dim2].m[i]) );
//			fname += ".dat";
//			save_dv_2d(fname, sobols[dim1], sobols[dim2]);
//
//			fname = "selected_dat_mse_Gauss/mse_" + to_string(sobols[dim1].d) + "_" + to_string(sobols[dim1].s) + "_" + to_string(sobols[dim1].a);
//			for (uint i = 0; i < sobols[dim1].s; ++i) fname += ( "_" + to_string(sobols[dim1].m[i]) );
//			fname += ("," + to_string(sobols[dim2].d) + "_" + to_string(sobols[dim2].s) + "_" + to_string(sobols[dim2].a) );
//			for (uint i = 0; i < sobols[dim2].s; ++i) fname += ( "_" + to_string(sobols[dim2].m[i]) );
//			fname += ".dat";
//			save_mse_2d(fname, mse_data_Gauss);
//		}
		return true;
	} else {
		return false;
	}
} // checkIntegrationGauss4D

bool checkStratificationAndIntegration4D(
		uint itrial,
        std::vector<SobolGenerator1D<sobolInt> >& sobols,
		const unsigned int dim1,
        const unsigned int dim2,
		const unsigned int dim3,
        const unsigned int dim4,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        std::ofstream& out,
        std::ofstream& out_strat,
		mt19937_64& gen,
        const unsigned int nTrials_integration,
		const bool owen_permut_flag,
		const bool max_tree_depth_32_flag,
        const bool dbg_flag) {

    bool test_strat4D = checkStratificationOne4D(sobols, dim1, dim2, dim3, dim4, octaveMin, octaveMax, dbg_flag);
    if (test_strat4D) {
        cout << " => "  << itrial << " : " << dim1 << " " << dim2 << " " << dim3 << " " << dim4 << " -> " << MSElimit_Gauss << " | " << kmax_best_found  << endl;

        std::uniform_int_distribution<uint32_t> unifFull(0);
        RNG randomGen;
        randomGen.seed(unifFull(gen));

        bool test_integration = checkIntegrationGauss4D(sobols, dim1, dim2, dim3, dim4, octaveMin, octaveMax, 1, nTrials_integration, randomGen,
    			owen_permut_flag, max_tree_depth_32_flag, dbg_flag);
    	if (test_integration) {
			out << sobols[dim1] << endl;
			out << sobols[dim2] << endl;
			out << sobols[dim3] << endl;
			out << sobols[dim4] << endl;
			out.flush();
			cout << " ############################################################################ " << endl;
			cout << sobols[dim1] << endl;
			cout << sobols[dim2] << endl;
			cout << sobols[dim3] << endl;
			cout << sobols[dim4] << endl;
			cout << endl << " ############################################################################ " << endl;
			cout.flush();
			return true;
		}
    }
    return false;
}	//checkStratificationAndIntegration4D

int main(int argc, char **argv) {
	srand48( time(NULL) );
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i) powers_of_two[i] = 1u << i;
    for(int i = 0; i < 16; ++i) powers_of_four[i] = 1u << (2*i);
	int dim1, dim2, dim3, dim4;
	unsigned int nTrials_integration = 32;
	std::string filename = "out.dat", filename_strat = "out_strat.dat";
    uint8_t NbThreadsMax = omp_get_max_threads();
	bool dbg_flag = false;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    bool all_test_are_OK;
	std::string dir_vectors_fname = "data/data_SobolPlusPlus/2D/upto_octave12/SobolPlusPlus_best_found.dat";
	unsigned int octaveMin = 2, octaveMax = 8;	// up to 16K spp only

	uint32_t seed = time(NULL); // 13374269;
	CLI::App app { "checkerStrat2D" };
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename );
//	app.add_option("--out_strat", filename_strat, "output filename_strat (ascii file), default: " + filename_strat );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: " + std::to_string(seed) );
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );
	app.add_option("--octaveMin", octaveMin, "min octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMin) );
	app.add_option("--octaveMax", octaveMax, "max octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMax) );
	app.add_option("--nTrials_integration", nTrials_integration, "nTrials_integration, default: " + std::to_string(nTrials_integration) );
	app.add_option("-p,--permut", owen_permut_flag, "owen_permut_flag, default: " + std::to_string(owen_permut_flag) );
	app.add_option("-m,--max_tree_depth_32", max_tree_depth_32_flag, "Owen or OwenPlus ? default: OwenPlus max_tree_depth_32_flag="+ std::to_string(max_tree_depth_32_flag)+"" );
	app.add_option("--MSElimit_Gauss", MSElimit_Gauss, "MSE goodness threshold during integration, default: " + std::to_string(MSElimit_Gauss) );
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + std::to_string(NbThreadsMax) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: " + std::to_string(dbg_flag) );
	CLI11_PARSE(app, argc, argv)

	srand48(seed);
	std::ofstream out(filename, std::ofstream::app);
	std::ofstream out_strat;
//	std::ofstream out_strat(filename_strat, std::ofstream::app);

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

    cout << "effective number of threads = "<< int(NbThreadsMax) << endl;
    cout << "octaveMin = "<< int(octaveMin) << "  " << "octaveMax = "<< int(octaveMax) << endl;
    cout << "MSElimit_Gauss = "<< MSElimit_Gauss << endl;

    cout << "Loading Sobol dir_vectors from " << dir_vectors_fname << endl;
    std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    loadSobolsFromFile(dir_vectors_fname, sobols, seed);		// read sobols from file
    cout << "Loading Sobol dir_vectors from " << dir_vectors_fname << " of size " << sobols.size() << endl;
    std::vector<SobolGenerator1D<uint32_t> > sobolsCpy = sobols;	// copy of sobols
    std::mt19937_64 gen( time(NULL) );

    int npairs = sobols.size() / 2;
    uint itrial = 0;
    while (true) {
    	itrial++;
        dim1 =  2 * floor(npairs*drand48());
        dim2 = dim1 + 1;
        dim3 = 2 * floor(npairs*drand48());
        dim4 = dim3 + 1;
        checkStratificationAndIntegration4D( itrial, sobols,
       		 dim1,dim2,dim3,dim4,octaveMin,octaveMax,
			 out,out_strat, gen,
			 nTrials_integration,owen_permut_flag,max_tree_depth_32_flag,dbg_flag);

	}
//
//    out.close();
////    out_strat.close();
//    if(all_test_are_OK) exit(0); else exit(1);
}
