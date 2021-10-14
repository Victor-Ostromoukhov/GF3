// optimizeScreenSpace2D.cpp -- vo Oct 2020
// based on solidangle2D.cpp (May 2018, Oct 2019)
// conventions: we call iptsXY 2D integer points (t_int_point2D) in the screen space (pixels)
//              we call ptsUV 2D floating points (t_point2D) in the sampling space (a priory, U and V are independent)
//				we call iseedsUV integer seeds (t_int_point2D) needed to generate ptsUV through our OwenPlus
// action:  1. read input_tile_seeds_128x128.dat and feel iseedsUV[i] of size iseedsUV containing seeds for dimensions dim1 and dim2
//          2. genrate ptsUV[i], using our OwenPlus, with owen_tree_depth=32
//          3. apply simulated annealing: at each step, swap one pair of seeds, and test whether the cost function decreases
//          4. write newly optimized seeds into output_tile_seeds_128x128.dat (we write each 10000 steps; expected N of iterations is above 10.000.000
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <list>
#include <math.h>       /* sqrt, pow */
#ifdef _MSC_VER
#include <process.h>
#else
#include <unistd.h>
#endif
#include <ctime>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <iterator>
#include <random>
#include <omp.h>

#include "../CLI11.hpp"
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"
#include "../includes/SeedTileHelper.hpp"

double drand48(void);
void srand48(long int seedval);

using namespace std;

//--------------------------------- constants

#define SIGMAI 2.1
#define SIGMAS 1.0
#define SIGMAISQ 4.41 //	(SIGMAI*SIGMAI)
#define SIGMASSQ 1 //(SIGMAS*SIGMAS)
#define NEIGHBORS_R 10		// 12: 436 neighbors, very safe;  11: 372 neighbors;  10: 305 neighbors

#define INIT_TEMPERATURE	.00001		//<<<<<<<<<< PAY ATTENTION HERE .00001 is probably the good value for 10.000.000 iterations

#define	TOTAL_N_CYCLES	10000000
#define OUTPUT_EACH_N	100000
#define DISPLAY_EACH_N	100

//--------------------------------- structures
struct t_int_point2D {
	unsigned long x, y;
};
struct t_point2D {
	float x, y;
};
struct t_point4D {
	float x, y, u, v;
};

//--------------------------------- global variable
int tile_size1d = 128;
int tile_size = 128 * 128;
#ifdef _MSC_VER
int pid = _getpid();
#else
int pid = getpid();
#endif
//--------------------------------- routines

float getEuclidDist(t_point2D pt1, t_point2D pt2) {
	float dx = pt2.x - pt1.x;
	float dy = pt2.y - pt1.y;
	return sqrt(dx * dx + dy * dy);
}	// getEuclidDist

int getIntToroidalDistSqWithWidth(t_int_point2D pt1, t_int_point2D pt2, int width) {
	long dx = abs((long) pt2.x - (long) pt1.x);
	long dy = abs((long) pt2.y - (long) pt1.y);
	if (dx > width / 2)
		dx = width - dx;
	if (dy > width / 2)
		dy = width - dy;
	return (dx * dx + dy * dy);
}	// getIntToroidalDistSqWithWidth

int getIntToroidalDistWithWidth(t_int_point2D pt1, t_int_point2D pt2, int width) {
	long dx = abs((long) pt2.x - (long) pt1.x);
	long dy = abs((long) pt2.y - (long) pt1.y);
	if (dx > width / 2)
		dx = width - dx;
	if (dy > width / 2)
		dy = width - dy;
	return sqrt(dx * dx + dy * dy);
}	// getIntToroidalDistWithWidth

void swap_mapping(int i, int *mapping) {
	unsigned int ind1, ind2, ind3, ind4, tmp, nswaps;
	ind1 = floor(tile_size * drand48());
	ind2 = floor(tile_size * drand48());
	tmp = mapping[ind1];
	mapping[ind1] = mapping[ind2];
	mapping[ind2] = tmp;
}	// swap_mapping

void write_tile(const string out_fname, const int *mapping, t_int_point2D *iptsXY, t_int_point2D *iseedsUV) {
	cout << "writing tile into " << out_fname << "\n";
	ofstream fp_out;
	fp_out.open(out_fname, ios::out);
	if (!fp_out.is_open()) {
		cout << out_fname << " cannot be written." << "\n";
		exit(1);
	};
	for (int i = 0; i < tile_size; i++) {
		t_int_point2D iptxy = iptsXY[i];
		t_int_point2D iseeduv = iseedsUV[mapping[i]];
		fp_out << iptxy.x << " " << iptxy.y << " " << iseeduv.x << " " << iseeduv.y << "\n";
	}
	fp_out.close();
}	// write_tile

void read_tile(const string in_fname, t_int_point2D *iptsXY, t_int_point2D *iseedsUV, const bool dbg_flag) {
	if (dbg_flag)
		cout << "reading tile from " << in_fname << endl;
	ifstream fp_in;
	fp_in.open(in_fname);
	if (!fp_in.is_open()) {
		cerr << in_fname << " cannot be read." << endl;
		exit(1);
	};
	int i = 0;
	do {
		if (i >= tile_size) {
			cerr << i << " something goes wrong in " << in_fname << endl;
			exit(1);
		}
		t_int_point2D iptxy, iseeduv;
		fp_in >> iptxy.x >> iptxy.y >> iseeduv.x >> iseeduv.y;
		iptsXY[i] = iptxy;
		iseedsUV[i] = iseeduv;
		i++;
//		cout << iptxy.x << " " << iptxy.y << " " << iseeduv.x << " " << iseeduv.y << endl;
	} while (!fp_in.eof());
	fp_in.close();
}	// write_tile

double get_1pt_Energy_nei(const int iRefPt, int *mapping, float **weightsxy, float **weightsuv, vector<int> *neighbors) {
	double totalNeiDist = 0.;
	for (int j = 0; j < neighbors[iRefPt].size(); j++) {
		totalNeiDist += weightsxy[iRefPt][neighbors[iRefPt][j]] * weightsuv[mapping[iRefPt]][mapping[neighbors[iRefPt][j]]];
	}
	return totalNeiDist;
}	// get_1pt_Energy_nei

double getTotalEnergy(int *mapping, float **weightsxy, float **weightsuv, vector<int> *neighbors) {
	double totalEnergy = 0.;
//#pragma omp parallel for reduction(+:totalEnergy)
#pragma omp parallel for schedule(dynamic) reduction(+:totalEnergy)
	for (int i = 0; i < tile_size; i++) { //tile_size
		totalEnergy += get_1pt_Energy_nei(i, mapping, weightsxy, weightsuv, neighbors);
	}
	return totalEnergy;
}

void savePPMfile(const std::string &filename, const int	size, const t_int_point2D *image) {
    std::ofstream ofs(filename,  std::ios::out);
    ofs << "P3" << endl;
    ofs << size << " " << size << endl;
    ofs << "255" << std::endl;
    for(int iy=0; iy < size; iy++) {
        for(int ix=0; ix < size; ix++) {
        	int ind = iy * size + ix;
            ofs << image[ind].x << " " << image[ind].y << " " << " 0 ";
        }
        ofs << endl;
    }
    ofs.close();
  }

string ToString(int value,int digitsCount) {
    ostringstream os;
    os<<setfill('0')<<setw(digitsCount)<<value;
    return os.str();
}

void loadDescDistance(std::vector<double> &descDistances,
                      const std::string filename,
                      const size_t size)
{
    double val;
    std::ifstream ifs(filename);
    for(auto idx=0; idx < size*size*size*size; ++idx)
    {
        ifs >> val;
        descDistances[idx] = val;
    }
    ifs.close();
}

int main(int argc, char **argv) {
	cout.precision(10);
	srand48((long int) time(NULL));
	double totalEnergy = 0;
	string in_fname = "";
	string out_fname = "data_ScreenSpace/tile_barycenters_128x128_ind1_2_iterXXXXXXXX.dat";
	uint8_t NbThreadsMax = omp_get_max_threads();
	int *current_mappingXY2UV;
	t_point2D *ptsUV;
	t_int_point2D *iptsXY, *iseedsUV, *iptsUV_res;
	float **weightsxy, **weightsuv;
	double temperature = INIT_TEMPERATURE;
	double kTermerature = 1. - 3. / 1000000.;
	unsigned int dim1 = 1, dim2 = -1;
	std::string dir_vectors_fname = "data/vo_sobol_init_tab.dat";
	uint32_t owen_tree_depth = OWENPLUS_TREE_DEPTH;
	bool dbg_flag = false;
    std::string importDescDistance_fname="";

	CLI::App app { "optimizeScreenSpace2D" };
	app.add_option("-i", in_fname, "input fname for tile_seeds");	//->required()->check(CLI::ExistingFile);
//	app.add_option("-o", out_fname, "output filename for tile_seeds");
	app.add_option("-t", temperature, "temperature (initial), default: " + to_string(temperature));
	app.add_option("-k", kTermerature, "temperature decreasing factor, default: " + to_string(kTermerature));
	app.add_option("-n,--nbThreadsMax", NbThreadsMax, "Maximum number of threads, default: max number of thread of the OS = " + to_string(NbThreadsMax));
	app.add_option("--dim1", dim1, "Sobol dim1 convetion: dim1=0 : van der Corput, default: " + to_string(dim1));
	app.add_option("--dim2", dim2, "Sobol dim12 default: default: " + to_string(dim2));
	app.add_option("-d,--dir_vectors", dir_vectors_fname, "dir_vectors_fname (ascii file), default: " + dir_vectors_fname);	// ->required()->check(CLI::ExistingFile);
	app.add_option("--dbg", dbg_flag, "dbg_flag, default: " + to_string(dbg_flag));
    app.add_option("--loadDescDistancesFromFile", importDescDistance_fname,"Import descriptor pairwise distance from file");
    CLI11_PARSE(app, argc, argv)
    if(dim2 == -1)  dim2 = dim1 + 1;
	if(in_fname == "") in_fname = "data_ScreenSpace/tile_barycenters_" + to_string(tile_size1d) + "x" + to_string(tile_size1d) +
							"_ind" + to_string(dim1) + "_" + to_string(dim2) + "_iter" + ToString(0, 8) + ".dat";

//	omp_set_dynamic(0);     // Explicitly disable dynamic teams
//	omp_set_num_threads(NbThreadsMax); // Use 32 threads for all consecutive parallel regions

	printf("================PROCESSING: %s size=%d temperature=%.10f\n", argv[0], tile_size1d, temperature);
	current_mappingXY2UV = (int*) malloc(sizeof(int) * tile_size);
	for (int i = 0; i < tile_size; i++) current_mappingXY2UV[i] = i;

// init of iptsXY & ptsUV
	cout << "Initialization..." << endl;
	ptsUV = (t_point2D*) malloc(sizeof(t_point2D) * tile_size);
	iptsXY = (t_int_point2D*) malloc(sizeof(t_int_point2D) * tile_size);
	iptsUV_res = (t_int_point2D*) malloc(sizeof(t_int_point2D) * tile_size);
	iseedsUV = (t_int_point2D*) malloc(sizeof(t_int_point2D) * tile_size);
//	neighbors = (std::vector<int> *) malloc(sizeof(vector<int> ) * tile_size);
	std::vector<int> *neighbors = new vector<int>[tile_size];

	SeedTile seedMap;
	seedMap.loadTile(in_fname);
	tile_size1d = seedMap.size;
	tile_size = tile_size1d * tile_size1d;

	if (seedMap.dimension != 2) {
		std::cerr << "[Error] Only optimizing 2D seeds for now, the seed map has dimension " << seedMap.dimension << std::endl;
		exit(2);
	}

	uint cpt = 0;
	for (auto i = 0; i < tile_size1d; ++i) {
		for (auto j = 0; j < tile_size1d; ++j) {
			t_int_point2D pix, uv;
			pix.x = i;
			pix.y = j;
			uv.x = seedMap(i, j, 0);
			uv.y = seedMap(i, j, 1);
			iptsXY[cpt] = pix;
			iseedsUV[cpt] = uv;
			++cpt;
		}
	}
	SeedTile tmpMap(seedMap.size, seedMap.dimension);

//	read_tile(in_fname, iptsXY, iseedsUV, dbg_flag);

	std::vector < SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
	loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file

	// init ptsUV
	cpt = 0;
	for (int i = 0; i < tile_size1d; ++i)
		for (int j = 0; j < tile_size1d; ++j) {
			t_point2D pts;
			pts.x = getOwenPlus1D_with_seed(sobols, dim1, 0, seedMap(i, j, 0));
			pts.y = getOwenPlus1D_with_seed(sobols, dim2, 0, seedMap(i, j, 1));
			ptsUV[cpt] = pts;
//			cout << i << " " << j << " | "  << seedMap(i, j, 0) << " "  << seedMap(i, j, 1) << " -> "  << pts.x << " " << pts.y << endl;
			++cpt;
		}

	// init of neighbors
	for (int i = 0; i < tile_size; i++) {
		neighbors[i].clear();
		t_int_point2D refptxy = iptsXY[i];
		for (int j = 0; j < tile_size; j++) {
//			if (i == j) continue;
			t_int_point2D ptxy = iptsXY[j];
			if ((getIntToroidalDistWithWidth(refptxy, ptxy, tile_size1d) < NEIGHBORS_R) && (1 != j)) {
				neighbors[i].push_back(j);
				if (dbg_flag)
					cout << i << "/" << j << " ref: " << refptxy.x << " " << refptxy.y << " | " << ptxy.x << " " << ptxy.y << " -> "
							<< getIntToroidalDistWithWidth(refptxy, ptxy, tile_size1d) << " -> " << neighbors[i].size() << endl;
			}
		}
	}
	if (dbg_flag)
		cout << "init neighbors: done." << endl;


    std::vector<double> descDistance( tile_size * tile_size,0.0);

    if (importDescDistance_fname!="")
        loadDescDistance(descDistance, importDescDistance_fname, tile_size1d);

	// init of weightsxy & weightsuv
	weightsxy = (float**) malloc(sizeof(float*) * tile_size);
	for (int i = 0; i < tile_size; i++)
		weightsxy[i] = (float*) malloc(sizeof(float) * tile_size);
	weightsuv = (float**) malloc(sizeof(float*) * tile_size);
	for (int i = 0; i < tile_size; i++)
		weightsuv[i] = (float*) malloc(sizeof(float) * tile_size);
	for (int i = 0; i < tile_size; i++) {
		t_int_point2D refptxy = iptsXY[i];
		t_point2D refptuv = ptsUV[i];
		float distxy, distuv;
		for (int j = 0; j < tile_size; j++) {
			t_int_point2D ptxy = iptsXY[j];
			distxy = (float) getIntToroidalDistSqWithWidth(refptxy, ptxy, tile_size1d);
			weightsxy[i][j] = exp(-distxy / SIGMAISQ);
			if (importDescDistance_fname == ""){
                t_point2D ptuv = ptsUV[j];
                distuv = getEuclidDist(refptuv, ptuv);
			} else {
			    distuv = descDistance[j + i * tile_size];
			}
			weightsuv[i][j] = exp(-distuv);
		}
	}
	if (dbg_flag) cout << "init weightsxy and weightsuv : done." << endl;
	cout << "Initialization: done." << endl;

	double initEnergy = getTotalEnergy(current_mappingXY2UV, weightsxy, weightsuv, neighbors);
	double bestFound = initEnergy;
	double lambda = kTermerature;
	for (int icycle = 0; icycle <= TOTAL_N_CYCLES; icycle++) {
		temperature = lambda * temperature;
		int swap_ind1 = floor(tile_size * drand48());
		int swap_ind2 = floor(tile_size * drand48());
		int tmp = current_mappingXY2UV[swap_ind1];
		current_mappingXY2UV[swap_ind1] = current_mappingXY2UV[swap_ind2];
		current_mappingXY2UV[swap_ind2] = tmp;
		totalEnergy = getTotalEnergy(current_mappingXY2UV, weightsxy, weightsuv, neighbors);
		double deltaE = (double) (totalEnergy - bestFound) / (double) totalEnergy;
		if (deltaE < 0) {
			bestFound = totalEnergy;
			if ((icycle % DISPLAY_EACH_N) == 0) {
				cout << "optimizeScreenSpace2D" << " " << icycle << "\t --- \t" << round(totalEnergy) << " " << round(bestFound) << "/" << round(initEnergy)
						<< " \t | " << round(1000000. * deltaE) << "\n";
			}
		} else {
			double expVal = exp(-deltaE / temperature);
			if (deltaE != 0 && drand48() < expVal) {
				if ((icycle % DISPLAY_EACH_N) == 0) {
					cout << "optimizeScreenSpace2D" << " " << icycle << "\t ^^^ \t" << round(totalEnergy) << " " << round(bestFound) << "/" << round(initEnergy)
							<< " \t | " << round(1000000. * deltaE) << " -> " << temperature << "\n";
					bestFound = totalEnergy;
				}
			} else {
				// swap_back
				if ((icycle % DISPLAY_EACH_N) == 0) {
					cout << "optimizeScreenSpace2D" << " " << icycle << "\t *** \t" << round(totalEnergy) << " " << round(bestFound) << "/" << round(initEnergy)
							<< " \t | " << "*" << " -> " << temperature << "\n";

				}
				int tmp = current_mappingXY2UV[swap_ind1];
				current_mappingXY2UV[swap_ind1] = current_mappingXY2UV[swap_ind2];
				current_mappingXY2UV[swap_ind2] = tmp;
			}
		}
		// output the result
		if ((icycle % OUTPUT_EACH_N) == 0) {
			if(icycle != 0) {
				for (auto cpt = 0; cpt < tile_size; ++cpt) {
					tmpMap(iptsXY[cpt].x, iptsXY[cpt].y, 0) = iseedsUV[current_mappingXY2UV[cpt]].x;
					tmpMap(iptsXY[cpt].x, iptsXY[cpt].y, 1) = iseedsUV[current_mappingXY2UV[cpt]].y;
				}
				out_fname = "data_ScreenSpace/tile_barycenters_" + to_string(tile_size1d) + "x" + to_string(tile_size1d) +
										"_ind" + to_string(dim1) + "_" + to_string(dim2) + "_iter" + ToString(icycle, 8) + ".dat";
				cout << "Saving tile into " << out_fname << endl;
				tmpMap.saveTile(out_fname);
			}
			for (auto i = 0; i < tile_size; i++) {
				iptsUV_res[i].x = floor(256. * ptsUV[current_mappingXY2UV[i]].x);
				iptsUV_res[i].y = floor(256. * ptsUV[current_mappingXY2UV[i]].y);
			}
			out_fname = "data_ScreenSpace/UVvaluses_" + to_string(tile_size1d) + "x" + to_string(tile_size1d) +
				"_ind" + to_string(dim1) + "_" + to_string(dim2) + "_iter" + ToString(icycle, 8) + ".ppm";
			cout << "Saving tile into " << out_fname << endl;
			savePPMfile(out_fname, tile_size1d, iptsUV_res);
			//std::cout<<"Saved energie= "<<getTotalEnergy(current_mappingXY2UV, weightsxy, weightsuv, neighbors)<<std::endl;
			//      write_tile(out_fname, current_mappingXY2UV, iptsXY, iseedsUV);
		}
	}
	return 0;
}

