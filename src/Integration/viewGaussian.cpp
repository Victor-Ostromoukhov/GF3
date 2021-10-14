//
// Created by lois on 23/11/2020.
//

#include <string>

#include "Integration/Integration.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "CLI11.hpp"

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
        cout << dir << endl;
        return 1.2 + std::sin(4 * M_PI * dir * point);
    }
    return 0;
}

int main(int argc, const char** argv){

    const int DIM = 2;
    string filename = "";
    string outfilename = "";
    int seed = 1234;
    int nbPts = 1;
    string dir_vectors_fname;
    vector<int> dims;
    int size = 1024;
    string funType = "gaussian";

    CLI::App app { "viewFunction" };
    app.add_option("-o", outfilename, "output filename for hdr image")->required();
    app.add_option("-g", seed, "Function index");
    app.add_option("--size", size, "Width of output image");
    app.add_option("--fun", funType, "Function type among: gaussian/slope/x/y/sine");
    CLI11_PARSE(app, argc, argv)

    std::vector<float> result(size*size);
    float mini = std::numeric_limits<float>::infinity();
//#pragma omp parallel for
    for(auto i=0; i < size; ++i) {
        for(auto j=0; j < size; ++j) {
            VecX<DIM> point(double(i) / size, double(j) / size);
            result[i + j * size] = getFunction(DIM, point, seed, funType);
            mini = std::min(result[i + j * size], mini);
        }
    }

    for (float &v : result){
        v -= mini;
    }

    stbi_write_hdr(outfilename.c_str(), size, size,1, result.data());

}