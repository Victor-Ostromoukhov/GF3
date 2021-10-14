// getGeneralizedL2Discrepancy_from_file.cpp

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <math.h>       /* fmod */
#include "CLI11.hpp"

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

struct type_pointND {double * coord ; };

//--------------------------------- global variables
double MSElimit = 3.;
uint32_t powers_of_two[CHAR_BIT * sizeof(sobolInt)];

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
	uint32_t seed = time(NULL);

	CLI::App app { "integrate2D_from_file" };
	app.add_option("-d,--nDims", nDims, "number of dimensions default: " + std::to_string(nDims) );
	app.add_option("-i,--input", input_filename, "input filename (ascii file), default: " + filename + ")" )->required();
	app.add_option("-o,--output", filename, "output filename (ascii file), default: " + filename + ")" );
	app.add_option("-f,--firstDim", firstDim, "d dimensions, starting from firstDim, default: "+ std::to_string(firstDim)+")" );
    app.add_option("-s,--seed", seed, "Random number generator seed. default: 13374269");
	app.add_option("--nbThreadsMax", NbThreadsMax, "Maximum number of threads (def: max number of thread of the OS = " + std::to_string(NbThreadsMax)+")");
	app.add_option("--integrandType", integrandType, "integrandType, possible values : 1 Heaviside 2 Gauss 3 Smooth 4 Cont 5 HeaviBell 6 HeaviCont 7 HeaviGauss default: " + std::to_string(integrandType) );
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: " + std::to_string(dbg_flag) );
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

    if(dbg_flag) {
    	for (int ipt = 0; ipt < points.size(); ipt++) {
    		for (unsigned int idim = 0; idim < nDims; idim++) {
    			cout << points[ipt].coord[idim] << " \t";
    		}
    		cout << endl;
    	}
    }

	double d = calculate_GeneralizedL2Discrepancy( nDims, points );
    cout << "calculate_GeneralizedL2Discrepancy for " << nDims <<"D " << points.size() << " pts -> " << d << endl;
    out << points.size() << " " <<d << endl;
    out.close();
}
