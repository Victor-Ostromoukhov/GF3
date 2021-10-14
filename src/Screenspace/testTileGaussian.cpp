//
// Created by lois on 23/11/2020.
//

#include <string>

#include "Integration/Integration.h"
#include "Samplers/SobolGenerator1D.h"
#include "Samplers/OwenScrambling.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "CLI11.hpp"
#include "includes/SeedTileHelper.hpp"

using namespace std;

template <typename Vec>
double getFunction(int DIM, const Vec& point, int seed, const string& type){
    if (type == "gaussian"){
        return getMultivariateGaussian(DIM, point.data(), tab_Gauss2D[seed % 16384].mu, tab_Gauss2D[seed % 16384].mxCInv);
    } else if (type == "x") {
        return point[0];
    } else if (type == "y") {
        return point[1];
    } else if (type == "slope") {
        mt19937_64 gen(seed);
        Vec dir = randomUnitVector<Vec>(DIM, gen);
        return 2 + point * dir;
    } else if (type == "sine") {
        mt19937_64 gen(seed);
        Vec dir = randomUnitVector<Vec>(DIM, gen);
        return 1.2 + std::sin(4 * M_PI * dir * point);
    }
    return 0;
}

int main(int argc, const char** argv){

    const int DIM = 2;
    string filename = "";
    string outfilename = "";
    int seed = 13374269;
    int nbPts = 1;
    string dir_vectors_fname;
    vector<int> dims;
    string funType = "gaussian";

    CLI::App app { "testTileGaussian" };
    app.add_option("-i", filename, "input fname for tile_seeds")->required()->check(CLI::ExistingFile);
    app.add_option("-o", outfilename, "output filename for hdr image")->required();
    app.add_option("-g", seed, "Gaussian index");
    app.add_option("-n", nbPts, "Number of point to use for integration");
    app.add_option("--fun", funType, "Function type among: gaussian/slope/x/y/sine");
    app.add_option("--dirs", dir_vectors_fname, "File name of the Sobol intialization table (e.g. ../../../data/sobol_init_tab.dat)")->required()->check(CLI::ExistingFile);
    app.add_option("--dims", dims, "Sobol dimensions to use")->required();

    CLI11_PARSE(app, argc, argv)

    SeedTile tile;
    tile.loadTile(filename);

    //init Sobol 01
    vector<SobolGenerator1D<uint32_t>> sobols;
    loadSobolsFromFile(dir_vectors_fname, sobols);

    std::cout << sobols[dims[0]] << endl;
    std::cout << sobols[dims[1]] << endl;

    //Init sobol points
    vector<array<uint32_t, DIM>> sobolPoint(nbPts);
    for (int i = 0; i < sobolPoint.size(); ++i){
        for (int j = 0; j < 2; ++j){
            sobolPoint[i][j] = sobols[dims[j]].getSobolInt(i);
        }
    }

    vector<VecX<DIM>> points(nbPts);
    int size = tile.size;
    std::vector<float> result(size*size);
//#pragma omp parallel for
    for(auto i=0; i < size; ++i) {
        for(auto j=0; j < size; ++j) {

            for (int indpt = 0; indpt < nbPts; ++indpt){
                for (int indDim = 0; indDim < DIM; ++indDim){
                    points[indpt][indDim] = double(OwenScrambling(sobolPoint[indpt][indDim], tile(i,j, indDim), 32)) / UINT32SOBOLNORM;
                }
            }
            double acc = 0.;
            for (const VecX<DIM>& point : points) {
                acc += getFunction(DIM, point, seed, funType);
            }
            acc /= nbPts;
            result[i+j*size] = acc;

        }
    }

    float mean = 0;
    float var = 0;
    for (const float &v : result){
        mean += v;
        var += v*v;
    }
    mean /= size*size;
    var /= size*size;
    var = abs(mean*mean - var);
    for (float &v : result){
        v -= mean;
        v *= 0.15 / sqrt(var);
        v += 0.3;
    }

    stbi_write_hdr(outfilename.c_str(), size, size,1, result.data());

}