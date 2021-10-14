// integrate2D_from_file.cpp

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

//--------------------------------- structures
typedef uint32_t sobolInt;

//--------------------------------- global variables
double MSElimit = 3.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];

struct type_pointND {double * coord ; };

int read_points_from_file(int nDims, istream& in, vector<type_pointND>& points){
	type_pointND pt_ND;
	pt_ND.coord = (double *)malloc(sizeof(double) * nDims);
	string line;
    points.clear();
    int count_in_pt_ND = 0, nlines = 0;;
    while(getline(in, line)){
        int c = line.find_first_not_of(" \t");
        if (line[c] != '#'){
            istringstream lineIn(line);
            double d = 0.;
            while (lineIn >> d){
            	pt_ND.coord[count_in_pt_ND] = d;
//                cout <<  count_in_pt_ND << " : " << d << " " << pt_ND.coord[count_in_pt_ND] << " | " << (count_in_pt_ND % nDims) << endl;
            	count_in_pt_ND++;
            	if (count_in_pt_ND % (nDims) == 0) {
            		points.push_back( pt_ND );
            		count_in_pt_ND = 0;
            		nlines++;
//                    cout << "-------------------- " << pt_ND.coord[0] << " " << pt_ND.coord[1] << " -> " << nlines << endl;
                    pt_ND.coord = (double *)malloc(sizeof(double) * nDims);
            	}
            }
         } else {
            break;
        }
    }
    return nlines;
}

//--------------------------------- routines related to integration
void integrate_One_PointSet(
		int firstDim,
		int nDims,
		vector<type_pointND>& in_pts,
		int npts,
        std::ofstream& out,
        const bool dbg_flag,
        const int integrandType = 1,
		const int nintegrands = 1024
		) {

    vector<double> integral_estimations;
    double mse_accumulator;

    vector<VecX<1>> points1(npts);
    vector<VecX<2>> points2(npts);
    vector<VecX<3>> points3(npts);
    vector<VecX<4>> points4(npts);
    vector<VecX<5>> points5(npts);
    vector<VecX<6>> points6(npts);
    vector<VecX<8>> points8(npts);
    vector<VecX<10>> points10(npts);
    vector<VecX<12>> points12(npts);

	mse_accumulator = 0.;
	for (int ipt = 0; ipt < npts; ipt++) {
		for (int idim = firstDim; idim < firstDim+nDims; idim++)
			switch(nDims) {
    		case 1:
    			points1[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 2:
    			points2[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 3:
    			points3[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 4:
    			points4[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 5:
    			points5[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 6:
    			points6[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 8:
    			points8[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 10:
    			points10[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
    		case 12:
    			points12[ipt][idim] = in_pts[ipt].coord[idim];
    			break;
			}
	}
	double this_mse;
	switch(nDims) {
	case 1:
		this_mse = calculate_mse(points1, integrandType, nintegrands);
		break;
	case 2:
		this_mse = calculate_mse(points2, integrandType, nintegrands);
		break;
	case 3:
		this_mse = calculate_mse(points3, integrandType, nintegrands);
		break;
	case 4:
		this_mse = calculate_mse(points4, integrandType, nintegrands);
		break;
	case 5:
		this_mse = calculate_mse(points5, integrandType, nintegrands);
		break;
	case 6:
		this_mse = calculate_mse(points6, integrandType, nintegrands);
		break;
	case 8:
		this_mse = calculate_mse(points8, integrandType, nintegrands);
		break;
	case 10:
		this_mse = calculate_mse(points10, integrandType, nintegrands);
		break;
	case 12:
		this_mse = calculate_mse(points12, integrandType, nintegrands);
		break;
	}
	double mean_mse = this_mse;
	if (dbg_flag) cout << "=====================> " << npts << " -> " << " " << mean_mse << endl;

	cout << std::setprecision(20) << npts << " \t" << mean_mse << endl;
	out  << std::setprecision(20) << npts << " \t" << mean_mse << endl;

} // integrate_One_PointSet

int main(int argc, char **argv) {
	srand48( time(NULL) );
    for(int i = 0; i < CHAR_BIT * sizeof(sobolInt); ++i){
        powers_of_two[i] = 1u << i;
    }


	unsigned int octaveMax = 14;

	std::string filename = "out.dat", input_filename = "in.dat";
	std::string dir_vectors_fname = "data/sobol_init_tab.dat";
    uint8_t NbThreadsMax = omp_get_max_threads();

	bool dbg_flag = false;
    bool owenPlus_permut_flag = true;
    bool checkIntegration_flag = true, checkStrat_flag = true;
    bool doubleCheck_flag = false;
    int integrandType = 1;
    unsigned int nDims=2, firstDim = 0;
	uint32_t seed = 13374269;
	int nintegrands = 1024;

	CLI::App app { "integrate2D_from_file" };
	app.add_option("-n,--nDims", nDims, "number of dimensions default: " + std::to_string(nDims) );
	app.add_option("-i,--input", input_filename, "input filename (ascii file), default: " + filename + ")" )->required();
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename + ")" );
	app.add_option("-f,--firstDim", firstDim, "d dimensions, starting from firstDim, default: "+ std::to_string(firstDim)+")" );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: 13374269");
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads (def: max number of thread of the OS = " + std::to_string(NbThreadsMax)+")");
	app.add_option("--integrandType", integrandType, "integrandType, possible values : 1 Heaviside 2 Gauss 3 Smooth 4 Cont 5 HeaviBell 6 HeaviCont 7 HeaviGauss default: " + std::to_string(integrandType) );
	app.add_option("--nintegrands", nintegrands, "number of integrands default: " + std::to_string(nintegrands) );
	app.add_option("-d,--dbg", dbg_flag, "dbg_flag, default: " + std::to_string(dbg_flag) );
	CLI11_PARSE(app, argc, argv)

	vector<type_pointND> points;
    std::ifstream in(input_filename);
	int npts_read = read_points_from_file(nDims, in, points);
    in.close();

	std::ofstream out(filename, std::ofstream::out);

    std::uniform_int_distribution<uint32_t> unifFull(0);
    std::mt19937_64 gen(seed );
    RNG randomGen;
    randomGen.seed(unifFull(gen));

    cout << "# integrateGauss2D_from_file integrandType = " << integrandType << "  for npts = " << npts_read << " output into " << filename << endl;

    integrate_One_PointSet(firstDim, nDims, points, npts_read, out, dbg_flag, integrandType, nintegrands);
    out.close();
}
