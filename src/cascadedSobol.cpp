#include <iostream>
#include <string>
#include <fstream>
#include <random>

#include "CLI11.hpp"

#include "Samplers/SobolGenerator1D.h"


using namespace std;

template< typename T >
inline void getUniformND( T *values,
                          const std::vector<SobolGenerator1D<uint32_t> >& sobols,
                          const uint32_t *seeds,
                          const int nDims,
                          const uint32_t n,
                          const uint32_t nbits,
                          const uint32_t owen_tree_depth = 32,
                          const bool owen_permut_flag = true )
{
    uint32_t indice = n;		// dim 0: take n-th point as index -> into 32-bit integer IDcode
    for(int idim = 0; idim < nDims; idim++)
    {
        indice = sobols[idim].getSobolInt(indice);	// radix-inversion + sobol

        uint32_t result = indice;
        if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
            result= OwenScrambling(result, seeds[idim], owen_tree_depth);
        values[idim]= T(result) / T(UINT32SOBOLNORM);	// final value (double) for this dimension

        indice = indice >> (32-nbits);				// this will be used as new index for the next dimension
    }
}

int main(int argc, char **argv) {


    CLI::App app { "cascadedSobol" };
    int nDims = 2;
    app.add_option("-d,--nDims", nDims, "number of dimensions to generate, default: "+ std::to_string(nDims)  );
    int npts = 1024;
    app.add_option("-n,--npts", npts, "number of points to generate, default: "+ std::to_string(npts) );
    int seed = 133742;
    app.add_option("-s,--seed", seed, "Random number generator seed. default: "+ std::to_string(seed)  );
    string input_dir_vectors = "data/vo_sobol_init_tab.dat";
    app.add_option("-i,--idv", input_dir_vectors, "input sobol initialisation (ascii file), default: " + input_dir_vectors );
    bool owen_permut_flag = false;
    app.add_option("-p,--permut", owen_permut_flag, "permut or not? when owen_permut_flag=0 ==> Sobol, default: "+ std::to_string(owen_permut_flag)  );
    bool dbg_flag = false;
    app.add_option("--dbg", dbg_flag, "dbg_flag, default: "+ std::to_string(dbg_flag)  );
    // Realizations and export "à la" utk
    int nbReal = 1;
    app.add_option("-m", nbReal, "number of realizations of the sampler, default: "+ std::to_string(nbReal));
    std::string output_fname = "out.dat";
    app.add_option("-o,--output", output_fname, "output fname, default: " + output_fname );
    CLI11_PARSE(app, argc, argv)

    //Create sobol matrices data
    std::vector<SobolGenerator1D<uint32_t> > sobols;	// array of sobol data per dim
    std::vector<uint32_t> d;
    std::vector<uint32_t> s;
    std::vector<uint32_t> a;
    std::vector<std::vector<uint32_t>> m;

    //Read sobol matrices value from file
    std::ifstream tableFile(input_dir_vectors);
    if (!tableFile.is_open()) {
        std::cerr << input_dir_vectors << " cannot be read." << std::endl;
        exit(1);
    };
    load_init_table(tableFile, d, s, a, m, nDims);
    if(dbg_flag){
        std::cout << "loading file " << input_dir_vectors << " " << d.size() << " entries." << std::endl;
    }
    init_sobols(sobols, d, s, a, m, seed, false);

    //Open output file
    std::ofstream ofs(output_fname , std::ios_base::out);

    //init random generator for owen scrambling
    std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);  // uniform distribution of uint32_ts between 0 and 2^32-1
    mt19937_64 gen(seed);

    auto nbits = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
    vector<double> points(npts * nDims);
    //For all realizations
    for(int realization=0; realization < nbReal; ++realization )
    {
        vector<uint32_t> realSeeds(nDims);
        for (int iDim = 0; iDim < nDims; ++iDim) {
            realSeeds[iDim] = unif32bits(gen);
        }

        for (int ipt = 0; ipt < npts; ipt++) {
            getUniformND(points.data() + ipt * nDims, sobols, realSeeds.data(), nDims, ipt, nbits, 32, owen_permut_flag);
        }

        //export points in ASCII
        for (int ipt = 0; ipt < npts; ipt++) {
            for (unsigned int idim = 0; idim < nDims; idim++) {
                ofs << std::setprecision(20) << std::fixed << points[ipt * nDims + idim] << " ";  // dims {0, ..., nDims-1}
            }
            ofs << endl;
        }
        //Diffent realizations are separated by #
        if (realization+1 != nbReal){
            ofs<<"#"<<std::endl;
        }
    }
    ofs.close();
}

