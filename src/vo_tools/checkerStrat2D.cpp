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
double MSElimit_Heaviside = 2., MSElimit_Gauss = 3.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];

bool checkStratificationOnePairOfDims(
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
		const unsigned int dim1,
		const unsigned int dim2,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
		const bool dbg_flag
		) {
	uint32_t x;
	uint32_t y;
	uint32_t shift;
	uint32_t npts1D;
	unsigned char squares[1024][1024];			// [2048][2048] for max octave = 22 [1024][1024] for max octave = 20, [128][128] for max octave = 14 [256][256] for max octave = 16

	for (uint iOctave = octaveMin; iOctave <= octaveMax; iOctave++) {	// almost useless
		uint npts = powers_of_two[iOctave];
		if ((iOctave & 1) == 0) { //enven octave
			npts1D = powers_of_two[(iOctave) / 2];
			shift = (N_SIGNIFICANT_BITS - (iOctave) / 2);
			// init squares
			for (uint iy = 0; iy < npts1D; iy++)
				for (uint ix = 0; ix < npts1D; ix++)
					squares[iy][ix] = 0;
			for (uint ipt = 0; ipt < npts; ipt++) {
				x = sobols[dim1].getSobolInt(ipt) >> shift;
				y = sobols[dim2].getSobolInt(ipt)  >> shift;
				if (squares[y][x] == 1) {
					return false;
				}
				squares[y][x] = 1;
			}
		} else { // odd octave
			npts1D = powers_of_two[(iOctave + 1) / 2];
			shift = (N_SIGNIFICANT_BITS - (iOctave + 1) / 2);
	        // init squares
	        for (uint iy = 0; iy < npts1D; iy++)
	            for (uint ix = 0; ix < npts1D; ix++)
	                squares[iy][ix] = 0;
			for (uint ipt = 0; ipt < npts; ipt++) {
				if (shift == N_SIGNIFICANT_BITS - 1) {
					x = sobols[dim1].getSobolInt(ipt) >> shift;
					y = 0;
				} else {
					x = sobols[dim1].getSobolInt(ipt) >> shift;
					y = sobols[dim2].getSobolInt(ipt) >> (shift + 1);
				}
				if (squares[y][x] == 1) {
					return false;
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
					y = sobols[dim2].getSobolInt(ipt) >> shift;
				} else {
					x = sobols[dim1].getSobolInt(ipt) >> (shift + 1);
					y = sobols[dim2].getSobolInt(ipt) >> shift;
				}
				if (squares[y][x] == 1) {
					return false;
				}
				squares[y][x] = 1;
			}
		}
	}
	return true;
}	// checkStratificationOnePairOfDims

bool checkIntegrationHeavisideAndGauss(
        const std::vector<SobolGenerator1D<sobolInt> >& sobols,
        const unsigned int dim1,
        const unsigned int dim2,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        const unsigned int octaveStep,
        const unsigned int nTrials_integration,
		RNG& gen,
		const bool checkHeaviside_flag,
        const bool doubleCheck_flag,
		const bool owen_permut_flag,
		const bool max_tree_depth_32_flag,
        const bool dbg_flag
		) {
	double ref_mse_Heaviside[21] = {0.2374293752557941061809, 0.04843594025943107755916,
			 0.01857461372233069896165, 0.0059587659848615330202,
			 0.002125090842474777926008, 0.0007425086509114782679433,
			 0.000262587504111541191875, 0.00009258146654830412688401,
			 0.00003274423002102671138814, 0.00001156717407949142021413,
			 4.100611832204145003691e-6, 1.446769681331099976097e-6,
			 5.117664729128487271652e-7, 1.809568557812224052382e-7,
			 6.390867465079974112358e-8, 2.26147770878309996805e-8,
			 7.997271132713943504418e-9, 2.83336900712258280052e-9,
			 9.964843461194076778478e-10, 3.528691374663081856636e-10,
			 1.25064682993350709765e-10};
	double ref_mse_Gauss[21] = {0.0796582128394425631468, 0.01694308125440395149108,
			 0.003196749500538342142475, 0.00054400118387911713079,
			 0.00008735387154042818054142, 0.00001356937905720574069797,
			 1.988368614877635147294e-6, 2.876193975477184105502e-7,
			 4.063750326475264807482e-8, 5.745095969331995386693e-9,
			 7.964993328250891814328e-10, 1.097177498432254959874e-10,
			 1.475214693978105045439e-11, 2.010715703683638913095e-12,
			 2.700196504339276144303e-13, 3.644218684478390299913e-14,
			 4.768256286709755607343e-15, 6.532619754070486065579e-16,
			 8.04796893326529949743e-17, 1.120110328333538985302e-17,
			 1.457758991951547084553e-18};
    vector<double> integral_estimations_Heaviside, integral_estimations_Gauss;
    vector<t_pointND> mse_data_Heaviside, mse_data_Gauss;
    double mse_accumulator_Heaviside, mse_accumulator_Gauss;
    unsigned int nTrials_integration_significant = nTrials_integration;;
	double kmax_Heaviside = 0., kmax_Gauss = 0.;
	if (dbg_flag)
#pragma omp critical (output)
		cout << "Entering integration test" << endl;
	for (int ioctave = octaveMin; ioctave <= octaveMax; ioctave += octaveStep) {
		size_t npts = powers_of_two[ioctave];
		uint32_t owen_tree_depth = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
		if(max_tree_depth_32_flag) owen_tree_depth = OWENPLUS_TREE_DEPTH;

		// init integral_estimations_Heaviside, integral_estimations_Gauss
		integral_estimations_Heaviside.clear();
		integral_estimations_Gauss.clear();
        vector<VecX<2>> points(npts);
		for (int iset = 0; iset < nTrials_integration; iset++) {
            uint32_t seed1 = gen();
            uint32_t seed2 = gen();
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
			}
			// calculate integral_estimations
			if(checkHeaviside_flag) {
				double this_mse_Heaviside = calculate_mse(points, INTEGRAND_TYPE_HEAVISIDE);		// in SobolPlusPlus/src/Integration/Integration.h
		        integral_estimations_Heaviside.push_back(this_mse_Heaviside);
			} else {
				integral_estimations_Heaviside.push_back(0.);
			}
			double this_mse_Gauss = calculate_mse(points, INTEGRAND_TYPE_GAUSS);		// in SobolPlusPlus/src/Integration/Integration.h
	        integral_estimations_Gauss.push_back(this_mse_Gauss);

		}
		mse_accumulator_Heaviside = 0.;
		mse_accumulator_Gauss = 0.;
	    for (int i = 0; i < nTrials_integration_significant ; i++) {
	    	if(checkHeaviside_flag) mse_accumulator_Heaviside += integral_estimations_Heaviside[i];
	    	mse_accumulator_Gauss += integral_estimations_Gauss[i];
	    }
	    if(dbg_flag) {
	    	for (int i = 0; i < nTrials_integration_significant ; i++) cout << "{" << integral_estimations_Heaviside[i] << " " << integral_estimations_Gauss[i] << "} " ;
	    	cout << endl;
	    }
		double mean_mse_Heaviside = mse_accumulator_Heaviside/(double)nTrials_integration_significant;
		double k_Heaviside = mean_mse_Heaviside / ref_mse_Heaviside[ioctave];
		kmax_Heaviside = max(kmax_Heaviside, k_Heaviside);
		double mean_mse_Gauss = mse_accumulator_Gauss/(double)nTrials_integration_significant;
		double k_Gauss = mean_mse_Gauss / ref_mse_Gauss[ioctave];
		kmax_Gauss = max(kmax_Gauss, k_Gauss);

		if (dbg_flag) {
			cout  << dim1 << "," << dim2 << " / " << sobols.size() << " | " << sobols[dim1].d - 1 << "," << sobols[dim2].d - 1 << " =================>>> " << points.size() << " pts <<<================= "
				 << nTrials_integration << " trials -> " << " " << mean_mse_Heaviside
				 << " \t k_Heaviside=" << k_Heaviside  << " \t k_Gauss=" << k_Gauss
				 << " \t kmax_Heaviside=" << kmax_Heaviside  << " \t kmax_Gauss=" << kmax_Gauss
				 << (kmax_Heaviside > MSElimit_Heaviside) << endl;
		}
		if (kmax_Heaviside > MSElimit_Heaviside || kmax_Gauss > MSElimit_Gauss) {
			return false;
		}
		if(doubleCheck_flag) {
			t_pointND this_measure;
			this_measure.coord[0] = npts;
			this_measure.coord[1] = mean_mse_Heaviside;
			mse_data_Heaviside.push_back( this_measure );
			this_measure.coord[1] = mean_mse_Gauss;
			mse_data_Gauss.push_back( this_measure );
		}
	}

	if (kmax_Heaviside < MSElimit_Heaviside && kmax_Gauss < MSElimit_Gauss) {
		if(doubleCheck_flag) {
			string fname;
			fname = "selected_dat_dv/dv_" + to_string(sobols[dim1].d) + "_" + to_string(sobols[dim1].s) + "_" + to_string(sobols[dim1].a);
			for (uint i = 0; i < sobols[dim1].s; ++i) fname += ( "_" + to_string(sobols[dim1].m[i]) );
			fname += ("," + to_string(sobols[dim2].d) + "_" + to_string(sobols[dim2].s) + "_" + to_string(sobols[dim2].a) );
			for (uint i = 0; i < sobols[dim2].s; ++i) fname += ( "_" + to_string(sobols[dim2].m[i]) );
			fname += ".dat";
			save_dv_2d(fname, sobols[dim1], sobols[dim2]);

			if(checkHeaviside_flag) {
				fname = "selected_dat_mse_Heaviside/mse_" + to_string(sobols[dim1].d) + "_" + to_string(sobols[dim1].s) + "_" + to_string(sobols[dim1].a);
				for (uint i = 0; i < sobols[dim1].s; ++i) fname += ( "_" + to_string(sobols[dim1].m[i]) );
				fname += ("," + to_string(sobols[dim2].d) + "_" + to_string(sobols[dim2].s) + "_" + to_string(sobols[dim2].a) );
				for (uint i = 0; i < sobols[dim2].s; ++i) fname += ( "_" + to_string(sobols[dim2].m[i]) );
				fname += ".dat";
				save_mse_2d(fname, mse_data_Heaviside);
			}

			fname = "selected_dat_mse_Gauss/mse_" + to_string(sobols[dim1].d) + "_" + to_string(sobols[dim1].s) + "_" + to_string(sobols[dim1].a);
			for (uint i = 0; i < sobols[dim1].s; ++i) fname += ( "_" + to_string(sobols[dim1].m[i]) );
			fname += ("," + to_string(sobols[dim2].d) + "_" + to_string(sobols[dim2].s) + "_" + to_string(sobols[dim2].a) );
			for (uint i = 0; i < sobols[dim2].s; ++i) fname += ( "_" + to_string(sobols[dim2].m[i]) );
			fname += ".dat";
			save_mse_2d(fname, mse_data_Gauss);
}
		return true;
	} else {
		return false;
	}
} // checkIntegrationHeavisideAndGauss

bool checkStratificationAndIntegration(
        std::vector<SobolGenerator1D<sobolInt> >& sobols,
		const unsigned int dim1,
        const unsigned int dim2,
        const unsigned int nTrialsdim1,
        const unsigned int nTrialsdim2,
        const unsigned int octaveMin,
        const unsigned int octaveMax,
        std::ofstream& out,
        std::ofstream& out_strat,
		mt19937_64& gen,
        const unsigned int nTrials_integration,
		const bool checkStrat_flag,
        const bool checkIntegration_flag,
        const bool checkHeaviside_flag,
        const bool doubleCheck_flag,
		const bool owen_permut_flag,
		const bool max_tree_depth_32_flag,
        const bool dbg_flag) {
    std::uniform_int_distribution<uint32_t> unif01(0, 1);
	bool test_checkStratification = true;
   	bool test_integration;
	for (unsigned int itrial1 = 1; itrial1 <= nTrialsdim1; itrial1++) {
#pragma omp critical (output)
        if (dbg_flag)
            cout << dim1 << "," << dim2 << "   \t" << itrial1 << "/" << nTrialsdim1 << endl;
        for (auto itrial2 = 1; itrial2 <= nTrialsdim2; itrial2++) {
        	test_checkStratification = true;
        	if(!doubleCheck_flag) {
				unsigned int second_entry = unif01(gen);
				sobols[dim1].m[1] = (second_entry << 1) ^ 1;
				sobols[dim2].m[1] = ((1 - second_entry) << 1) ^ 1;
				for (int k = 2; k < sobols[dim1].s; ++k) {
					//Random in [0, 2^k[
					std::uniform_int_distribution<uint32_t> unifFull(0, powers_of_two[k] - 1);
					unsigned int i = unifFull(gen);
					sobols[dim1].m[k] = (i << 1) ^ 1;
				}
				for (int k = 2; k < sobols[dim2].s; ++k) {
					std::uniform_int_distribution<uint32_t> unifFull(0, powers_of_two[k] - 1);
					unsigned int i = unifFull(gen);
					sobols[dim2].m[k] = (i << 1) ^ 1;
				}
				sobols[dim1].generateMatrix();
				sobols[dim2].generateMatrix();
				if (dbg_flag)
#pragma omp critical (output)
					{
					cout << "=============== checking for\n";
                    cout << sobols[dim1] << endl;
                    cout << sobols[dim2] << endl;
					}
        	} else {
				if (dbg_flag) {
#pragma omp critical (output)
					{
					cout << "=============== double checking for\n";
                    cout << sobols[dim1] << endl;
                    cout << sobols[dim2] << endl;
					}
				}
			}
        	if (checkStrat_flag) test_checkStratification = checkStratificationOnePairOfDims(sobols, dim1, dim2, octaveMin, octaveMax, dbg_flag);
        	if(checkStrat_flag && !test_checkStratification) continue;
            if (test_checkStratification) {
                //We have found a good pair of dir vectors
#pragma omp critical (output)
                {
                    cout << dim1 << "," << dim2 << "   \t" << itrial1 << "/" << nTrialsdim1 << "   \t" << itrial2 << "/" << nTrialsdim2 << endl;
//                    cout << sobols[dim1] << endl;
//                    cout << sobols[dim2] << endl;
                    cout.flush();
//                    out_strat << sobols[dim1] << endl;
//                    out_strat << sobols[dim2] << endl;
//                    out_strat.flush();
                }
                std::uniform_int_distribution<uint32_t> unifFull(0);
                RNG randomGen;
                randomGen.seed(unifFull(gen));
                if (checkIntegration_flag)
                	test_integration = checkIntegrationHeavisideAndGauss(sobols, dim1, dim2, octaveMin, octaveMax, 1, nTrials_integration, randomGen,
                			checkHeaviside_flag, doubleCheck_flag, owen_permut_flag, max_tree_depth_32_flag, dbg_flag);
                if (checkIntegration_flag && test_integration) {    // check with octaveStep=1
#pragma omp critical (outputFile)
                    {
                        out << sobols[dim1] << endl;
                        out << sobols[dim2] << endl;
                        out.flush();
                        cout << " ############################################################################ " << endl;
                        cout << sobols[dim1] << endl;
                        cout << sobols[dim2] << endl;
                        cout << endl << " ############################################################################ " << endl;
                        cout.flush();
//                        if (dim2 == dim1+1) exit(0);
                    }
                }
            }
        }
    }
    return (test_checkStratification && test_integration);
} // checkStratificationAndIntegration

int main(int argc, char **argv) {
	srand48( time(NULL) );
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i) powers_of_two[i] = 1u << i;
	int dim1=1;
    int dim2=-1;
    unsigned int nTrialsdim1 = 64*1024;
    unsigned int nTrialsdim2 = 1024; // 4096 8192 16384 32768 65536
	unsigned int octaveMin = 2, octaveMax = 8;	// up to 16K spp only
	unsigned int nTrials_integration = 32;
	std::string filename = "out.dat", filename_strat = "out_strat.dat";
	std::string dir_vectors_fname = "data/sobol_init_tab.dat";
    uint8_t NbThreadsMax = omp_get_max_threads();
	bool dbg_flag = false;
    bool owen_permut_flag = true, max_tree_depth_32_flag = true;
    bool checkIntegration_flag = true, checkStrat_flag = true, checkHeaviside_flag = false;	// checkHeaviside_flag only when double-checking
    bool doubleCheck_flag = false;
    bool all_test_are_OK;

	uint32_t seed = time(NULL); //13374269;
	CLI::App app { "checkerStrat2D" };
	app.add_option("-i,--dim1", dim1, "Index of the first Sobol pair (==dim1) convention: dim1=0 : van der Corput")->required();
	app.add_option("-j,--dim2", dim2, "Index of the second Sobol pair (==dim2) default: dim2=dim1+1") ;
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename );
//	app.add_option("--out_strat", filename_strat, "output filename_strat (ascii file), default: " + filename_strat );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: " + std::to_string(seed) );
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname );
	app.add_option("--octaveMin", octaveMin, "min octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMin) );
	app.add_option("--octaveMax", octaveMax, "max octave value (will check all 2^i spp for i {octaveMin..octaveMax}), default: " + std::to_string(octaveMax) );
	app.add_option("--nTrials_dv_dim1", nTrialsdim1, "nTrials_dv_dim1, default: " + std::to_string(nTrialsdim1) );
	app.add_option("--nTrials_dv_dim2", nTrialsdim2, "nTrials_dv_dim2, default: " + std::to_string(nTrialsdim2) );
	app.add_option("--nTrials_integration", nTrials_integration, "nTrials_integration, default: " + std::to_string(nTrials_integration) );
	app.add_option("-p,--permut", owen_permut_flag, "owen_permut_flag, default: " + std::to_string(owen_permut_flag) );
	app.add_option("-m,--max_tree_depth_32", max_tree_depth_32_flag, "Owen or OwenPlus ? default: OwenPlus max_tree_depth_32_flag="+ std::to_string(max_tree_depth_32_flag)+"" );
	app.add_option("--checkStrat", checkStrat_flag, "checkStrat, default: " + std::to_string(checkStrat_flag) );
	app.add_option("--checkIntegration", checkIntegration_flag, "checkIntegration_flag, default: " + std::to_string(checkIntegration_flag) );
	app.add_option("--checkHeaviside", checkHeaviside_flag, "checkHeaviside_flag, default: " + std::to_string(checkHeaviside_flag) );
	app.add_option("--MSElimit_Heaviside", MSElimit_Heaviside, "MSE goodness threshold during integration, default: " + std::to_string(MSElimit_Heaviside) );
	app.add_option("--MSElimit_Gauss", MSElimit_Gauss, "MSE goodness threshold during integration, default: " + std::to_string(MSElimit_Gauss) );
	app.add_option("--doubleCheck", doubleCheck_flag, "doubleCheck_flag (use direction vectors from file + output mse data), default: " + std::to_string(doubleCheck_flag) );
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + std::to_string(NbThreadsMax) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: " + std::to_string(dbg_flag) );
	CLI11_PARSE(app, argc, argv)

	std::ofstream out(filename, std::ofstream::app);
	std::ofstream out_strat;
//	std::ofstream out_strat(filename_strat, std::ofstream::app);

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

    if (dim2 != -1){
        mt19937_64 gen(seed);
        if(nTrialsdim1 == 1 && nTrialsdim2 == 1) {
            std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
            if(dbg_flag) cout << "Loading Sobol dir_vectors from " << dir_vectors_fname << endl;
            loadSobolsFromFile(dir_vectors_fname, sobols, seed);		// read sobols from file
            all_test_are_OK = checkStratificationAndIntegration(sobols,
                                dim1,
                                dim2,
                                nTrialsdim1,
                                nTrialsdim2,
                                octaveMin,
                                octaveMax,
                                out,
                                out_strat,
                                gen,
                                nTrials_integration,
    							checkStrat_flag,
                                checkIntegration_flag,checkHeaviside_flag,
                                doubleCheck_flag,
								owen_permut_flag,
								max_tree_depth_32_flag,
                                dbg_flag);
        } else {
#pragma omp parallel for
        	for (int trial = 0; trial <= 64; ++trial) {
        	    std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
        	    if(dbg_flag) cout << "Loading Sobol dir_vectors from " << dir_vectors_fname << endl;
        	    loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file
            	all_test_are_OK = checkStratificationAndIntegration(sobols,
                                dim1,
                                dim2,
                                nTrialsdim1,
                                nTrialsdim2,
                                octaveMin,
                                octaveMax,
                                out,
                                out_strat,
                                gen,
                                nTrials_integration,
        						checkStrat_flag,
                                checkIntegration_flag,checkHeaviside_flag,
                                doubleCheck_flag,
        						owen_permut_flag,
        						max_tree_depth_32_flag,
                                dbg_flag);
        	}
        }

        if(all_test_are_OK) exit(0); else exit(1);
    }

    int dim2from = dim1 + 1;
	int dim2to = max(dim1 + 64, 200);
	if (dim1 == 0) dim2from = 2;
//	if (dim1 <= 7) dim2from = 19;
//	if (dim1 <= 12) dim2from = 17;
//	if (dim1 <= 9) dim2from = 11;
//	if (dim1 <= 7) dim2from = 19;
//	if (dim1 <= 4) dim2from = 17;
//	if (dim1 <= 2) dim2from = 21;

    cout << "MSElimit_Gauss = " << MSElimit_Gauss << endl;
    cout << "effective number of threads = " << int(NbThreadsMax) << endl;
//#pragma omp parallel for schedule(dynamic)
#pragma omp parallel for
	for (int dim2 = dim2from; dim2 <= dim2to; ++dim2) {
	    std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
	    if(dbg_flag) cout << "Loading Sobol dir_vectors from " << dir_vectors_fname << endl;
	    loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file

	    std::vector<SobolGenerator1D<uint32_t> > sobolsCpy = sobols;	// copy of sobols
        std::mt19937_64 gen(seed + 76543*dim1 + 54321*dim2);

#pragma omp critical (output)
        cout << "Checking Sobol indices (" << dim1 << "," << dim2 << ") for nTrialsdim1="
        		<< nTrialsdim1 << " nTrialsdim2=" << nTrialsdim2 << " and npts between 2^" << octaveMin << " and 2^" << octaveMax << " output into " << filename
				<< endl;
//#pragma omp parallel
        all_test_are_OK = checkStratificationAndIntegration(sobolsCpy,
		                    dim1,
		                    dim2,
		                    nTrialsdim1,
		                    nTrialsdim2,
                            octaveMin,
                            octaveMax,
                            out,
                            out_strat,
		                    gen,
		                    nTrials_integration,
							checkStrat_flag,
                            checkIntegration_flag,checkHeaviside_flag,
		                    doubleCheck_flag,
							owen_permut_flag,
							max_tree_depth_32_flag,
		                    dbg_flag);

	}

    out.close();
//    out_strat.close();
    if(all_test_are_OK) exit(0); else exit(1);
}
