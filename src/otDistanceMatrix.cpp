//
// Created by lois on 13/11/2020.
//

#include <vector>
#include <array>
#include <string>
#include <random>
#include <fstream>
#include <chrono>
#include "Transport/SlicedOptimalTransportNBall.h"
#include "Samplers/OwenScrambling.h"
#include "Samplers/SobolGenerator1D.h"
#include "includes/SeedTileHelper.hpp"
#include "Transport/network_simplex_simple.h"
#include "CLI11.hpp"

using namespace std;

#define DIM 2

typedef int arc_id_type; // {short, int, int64_t} ; Should be able to handle (n1*n2+n1+n2) with n1 and n2 the number of nodes (INT_MAX = 46340^2, I64_MAX = 3037000500^2)
typedef int supply_type; // {float, double, int, int64_t} ; Should be able to handle the sum of supplies and *should be signed* (a demand is a negative supply)
typedef double cost_type;  // {float, double, int, int64_t} ; Should be able to handle (number of arcs * maximum cost) and *should be signed*

struct TsFlow {
    int from, to;
    supply_type amount;
};

template<typename T>
cost_type sqrdistance2d(T* a, T* b) {
    return (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]);
}

double computeOTCost(const vector<VecX<DIM>> & points1, const vector<VecX<DIM>> &points2){
    using namespace lemon;

    typedef FullBipartiteDigraph Digraph;
    DIGRAPH_TYPEDEFS(FullBipartiteDigraph);

    arc_id_type n1 = points1.size(), n2 = points2.size(); // just demo for the case n1=n2 ; adapt otherwise
    std::vector<supply_type> weights1(n1, 1), weights2(n2, -1); // works faster with integer weights though

    Digraph di(n1, n2);
    NetworkSimplexSimple<Digraph, supply_type, cost_type, arc_id_type> net(di, true, n1 + n2, n1*n2);

    //Set costs
    arc_id_type idarc = 0;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            cost_type d = sqrdistance2d(points1[i].data(), points2[j].data());
            net.setCost(di.arcFromId(idarc), d);
            idarc++;
        }
    }

    //Set weights
    net.supplyMap(weights1.data(), n1, weights2.data(), n2);

    net.run();

    return net.totalCost();

}

void addNextPoints(vector<VecX<DIM>> &points, const vector<array<uint32_t, DIM>> &sobolPoint, int nbPts, const SeedTile &tile, int i, int j){

    size_t oldSize = points.size();
    points.resize(nbPts);
    for (int indpt = oldSize; indpt < nbPts; ++indpt){
        for (int indDim = 0; indDim < DIM; ++indDim){
            points[indpt][indDim] = double(OwenScrambling(sobolPoint[indpt][indDim], tile(i,j, indDim), 32)) / UINT32SOBOLNORM;
        }
    }
}

int main(int argc, const char **argv){

    srand(0);

    string outputFile = "out.dat";
    int size = 128;
    int seed = 13374269;
    int nbOctaves = 6;
    int startOctave = 0;

    CLI::App app { "otCostCompare" };
    app.add_option("-o", outputFile, "Output file. Default:out.dat");
    app.add_option("--start", startOctave, "Start octave default 0 (1 spp)");
    app.add_option("--end", nbOctaves, "End octave default 6 (64 spp)");
    app.add_option("-n", size, "Size of generated square tile. Default:128");
    app.add_option("-s", seed, "Seed for random generator. Default:13374269");
    CLI11_PARSE(app, argc, argv)

	    
    std::cerr << "Initialisation..." << std::endl;

    //init Sobol 01
    vector<SobolGenerator1D<uint32_t>> sobol(2);
    sobol[0] = SobolGenerator1D<uint32_t>(0, 0, 0, std::vector<uint32_t>(1));
    sobol[1] = SobolGenerator1D<uint32_t>(2, 1, 0, std::vector<uint32_t>(1, 1));


    //Init sobol points
    vector<array<uint32_t, DIM>> sobolPoint(1u << nbOctaves);
    for (int i = 0; i < sobolPoint.size(); ++i){
        for (int j = 0; j < 2; ++j){
            sobolPoint[i][j] = sobol[j].getSobolInt(i);
        }
    }

    //Allocate Points and distance tabs
    int nbInits = size * size;
    vector<vector<VecX<DIM>>> points(nbInits);
    vector<vector<double>> distances(nbOctaves - startOctave + 1);
    for (vector<double>& tab : distances){
        tab.resize(nbInits * nbInits);
    }

    //Init seeds
    SeedTile tile(size, DIM);
    mt19937_64 gen(seed);
    uniform_int_distribution<uint32_t> unif(0);
    for (uint32_t &v : tile.seedmap){
        v = unif(gen);
    }

    std::cerr << "Calcul des distances..." << std::endl;

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    for (int indOct = startOctave; indOct <= nbOctaves; ++indOct){
	std::cerr << "\tOctave numero: " << indOct << std::endl; 
        unsigned int nbPts = 1u << indOct;
        for (int indY = 0; indY < size; ++indY){
            for (int indX = 0; indX < size; ++indX){
                addNextPoints(points[indX + indY * size], sobolPoint, nbPts, tile, indX, indY);
            }
        }
#pragma omp parallel for schedule(dynamic)
        for (int indSeed1 = 0; indSeed1 < nbInits; ++indSeed1){
            for (int indSeed2 = indSeed1 + 1; indSeed2 < nbInits; ++indSeed2){
                //double dist = computeDiscreteSlicedOTCost(points[indSeed1], points[indSeed2], 2048, 1, gen);
                double dist = computeOTCost(points[indSeed1], points[indSeed2]);
                distances[indOct - startOctave][indSeed2 + indSeed1 * nbInits] = dist;
                distances[indOct - startOctave][indSeed1 + indSeed2 * nbInits] = dist;
            }
        }
        double mean = 0.;
        double var = 0.;
        int nbValues = 0;
        for (int indSeed1 = 0; indSeed1 < nbInits; ++indSeed1){
            for (int indSeed2 = indSeed1 + 1; indSeed2 < nbInits; ++indSeed2){
                nbValues += 1;
                double value = distances[indOct - startOctave][indSeed1 + indSeed2 * nbInits];
                double dvalue = (value - mean) / nbValues;
                mean += dvalue;
                var += dvalue*(value - mean);
            }
        }
        var /= nbValues-1;
        double eq = std::sqrt(var);
        for (double &v : distances[indOct - startOctave]){
            v -= mean;
            v /= eq;
        }
    }

    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    std::cerr<<"time: "<< (double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / 1000000.<<std::endl;
    std::cerr << "Ecriture des fichiers" << std::endl;
//    cout << "seeds:" << endl;
    tile.saveTile(outputFile);

    vector<double> total(nbInits * nbInits);
//    cout << "distances:" << endl;
    for (int indOctave = startOctave; indOctave <= nbOctaves; ++indOctave){
//        cout << "octave " << indOctave << ":" << endl;
        for (int indSeed = 0; indSeed < nbInits; ++indSeed) {
            for (int indSeed2 = 0; indSeed2 < nbInits; ++indSeed2) {
                total[indSeed + indSeed2 * nbInits] += distances[indOctave - startOctave][indSeed + indSeed2 * nbInits];

//                cout << distances[indOctave][indSeed + indSeed2 * seeds.size()] << " ";
            }
//            cout << endl;
        }
//        cout << endl;
    }
    double val = total[0];
    for (double &v : total){
        v -= val;
    }

    //cout << "total:" << endl;
    for (int indSeed = 0; indSeed < nbInits; ++indSeed) {
        for (int indSeed2 = 0; indSeed2 < nbInits; ++indSeed2) {
            cout << total[indSeed + indSeed2 * nbInits] << " ";
        }
        cout << endl;
    }


}
