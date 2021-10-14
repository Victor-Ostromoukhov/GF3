#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>

#include <CLI11.hpp>
#include "../Samplers/OwenScrambling.h"
#include "../Samplers/SobolGenerator1D.h"

using namespace std;

#define NPIXELS 256
#define NSPP 1024

struct t_int_seed2D {
	uint32_t dim1, dim2;
};

struct t_int_point2D {
	uint32_t x, y;
};

struct t_point2D {
	float u, v;
};

inline float get_pixel_sample(
		const std::vector<SobolGenerator1D<uint32_t> >& sobols,
		const t_int_point2D pixel,
		const uint32_t spp,
		const uint32_t dim,
		const uint32_t owen_tree_seed
		) {
	return getOwenPlus1D_with_seed(sobols, dim, spp, owen_tree_seed); 	// getOwenPlus1D_with_seed() is in src/Samplers/SobolGenerator1D.h
}

std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);	// uniform distribution of integers between 0 and 2^32-1
std::uniform_int_distribution<uint32_t> unif_pixels(0, NPIXELS-1);
std::uniform_int_distribution<uint32_t> unif_spp(1, NSPP);
std::uniform_int_distribution<uint32_t> unif_01(0, 1);

int main(int argc, char **argv)
{
  CLI::App app { "demo owen plus 2D" };
	unsigned int dim1 = 3, dim2 = 5, npts = 16;
  std::string dir_vectors_fname; // = "../../../data/sobol_init_tab.dat";
  app.add_option("-d,--dirs", dir_vectors_fname, "File name of the Sobol intialization table (e.g. ../../../data/sobol_init_tab.dat)")->required()->check(CLI::ExistingFile);
  CLI11_PARSE(app, argc, argv)

    std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    loadSobolsFromFile(dir_vectors_fname, sobols);		// read sobols from file and fill appropriate structures

    // first, output 16 2-d points (dim1 & dim2)
	for (int ipt = 0; ipt < npts; ipt++) {
		float x = getOwenPlus1D(sobols, dim1, ipt);	// getOwenPlus1D() is in src/Samplers/SobolGenerator1D.h
		float y = getOwenPlus1D(sobols, dim2, ipt);
		cout << "pt = " << ipt << " : \t" << std::setprecision(16) << x << "," << y << endl;
	}

    // now, output random-access per-pixel sample values, up to nspp
	t_int_seed2D per_pixel_seeds[NPIXELS][NPIXELS];
	// init per_pixel_seeds[][]
	uint32_t seed = 13374269;
	std::mt19937_64 gen(seed);
	for (int ipix1 = 0; ipix1 < NPIXELS; ipix1++) {
		for (int ipix2 = 0; ipix2 < NPIXELS; ipix2++) {
			per_pixel_seeds[ipix2][ipix1].dim1 = unif32bits(gen);
			per_pixel_seeds[ipix2][ipix1].dim2 = unif32bits(gen);
		}
	}
	// random-access to pixels/spp/dim
	for (int itrial = 0; itrial < 10; itrial++) {
		t_int_point2D pixel = {unif_pixels(gen), unif_pixels(gen)};
		uint32_t spp = unif_spp(gen);
		uint32_t dim = dim1;
		uint32_t seed = per_pixel_seeds[pixel.y][pixel.x].dim1;
		if(unif_01(gen) == 0) {
			dim = dim2;
			seed = per_pixel_seeds[pixel.y][pixel.x].dim2;
		}
		float value = get_pixel_sample(sobols, pixel, spp, dim, seed);
		cout << itrial << " pixel=(" << pixel.x << "," << pixel.y << ") spp=" << spp << " dim=" << dim << " -> " << value << endl;
	}

}
