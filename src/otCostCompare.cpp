//
// Created by lpaulin on 17/09/20.
//
#include "CLI11.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#ifndef _MSC_VER
#include <sys/time.h>
#endif
#include <random>

#include <geogram/basic/common.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/numerics/predicates.h>

#include "Samplers/OwenScrambling.h"
#include "Samplers/SobolGenerator1D.h"
#include "Points/VecX.h"
#include "Points/mapping.h"
#include "Transport/SlicedOptimalTransportNBall.h"
#include "Transport/Wasserstein.h"
#include "Integration/Gaussians.h"

template <int N>
std::vector<double> testSobol(double ( *fun) (const std::vector<VecX<N>>&), const std::vector<SobolGenerator1D<uint32_t>> & sobol, uint32_t nbPts, uint32_t nbReals, uint32_t seed){

    RNG rand;
    rand.seed(seed);
    std::vector<double> accumulator(nbReals);

     //Shuffle nbReals times and average fun
    for (int idReal = 0; idReal < nbReals; ++idReal){

        uint32_t localSeed = rand();
        //Generate Sobol points
        std::vector<VecX<N>> pts(nbPts);
        for(uint32_t idPts = 0; idPts < nbPts; ++idPts){
            for (uint32_t idDim = 0; idDim < N; ++idDim){
                pts[idPts][idDim] = OwenScrambling(sobol[idDim].getSobolInt(idPts), localSeed, OWENPLUS_TREE_DEPTH) / double(UINT32SOBOLNORM);
                //pts[idPts][idDim] = sobol[idDim].getSobolDouble(idPts);
            }
std::cout << std::endl;
        }

        accumulator[idReal] = fun(pts);

    }

    return accumulator;

}


template <int N>
std::vector<double> testWN(double ( *fun) (const std::vector<VecX<N>>&, uint32_t), uint32_t nbPts, uint32_t nbReals, uint32_t seed){

    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> unif(0.,1.);
    std::vector<double> accumulator(nbReals);

    //Shuffle nbReals times and average fun
    for (int idReal = 0; idReal < nbReals; ++idReal){

        //Generate WN points
        std::vector<VecX<N>> pts(nbPts);
        for(uint32_t idPts = 0; idPts < nbPts; ++idPts){
            for (uint32_t idDim = 0; idDim < N; ++idDim){
                pts[idPts][idDim] = unif(gen);
            }
        }
        accumulator[idReal] = fun(pts, seed);
    }

    return accumulator;

}


template <int N>
double computeSOTCost(const std::vector<VecX<N>>& pts, uint32_t seed){
    static const int nb_threads = 1;
    static const int nbSlices = 2048;
    static std::mt19937_64 gen(seed);

    std::vector<VecX<N> > mappedPoints(pts);
    for (VecX<N>& v : mappedPoints){
        v = NCube2NBall(2 * v - 1);
    }
    return computeSlicedOTCost(mappedPoints, nbSlices, nb_threads, gen);
}

template <int N>
double computeIntegration(const std::vector<VecX<N>>& pts, uint32_t seed){
    return calculate_mse_Gaussian(pts, 5, 1024, seed);
}


int main(int argc, const char** argv){

    std::string inputFile = "../../data/sobol_init_tab.dat"; //"/home/lpaulin/Documents/Projects/SobolPlusPlus/data/data_SobolPlusPlus_strat/SobolPlusPlus_best_found_strat.dat"; //
    uint32_t seed = 13374269;
    uint32_t nbPairs = 1;
    uint32_t offset = 0;
    uint32_t octaveMax = 12;
    uint32_t octaveDepart = 2;
    uint32_t nbptsDepart = 4;
    uint32_t nbOwens = 1;

    //Gestion des parametres
    CLI::App app { "otCostCompare" };
    app.add_option("-i", inputFile, "Init vectors file");
    app.add_option("--offset", offset, "Offset in input file");
    app.add_option("-n", nbPairs, "Number of pairs to process");
    app.add_option("-s", seed, "Seed for random generator. Default:13374269");
    CLI11_PARSE(app, argc, argv)

    //Initialisation de la table Sobol
    std::vector<uint32_t> d(1,0);
    std::vector<uint32_t> s(1,0);
    std::vector<uint32_t> a(1,0);
    std::vector<std::vector<uint32_t>> m(1);
    /*
    std::vector<uint32_t> d;
    std::vector<uint32_t> s;
    std::vector<uint32_t> a;
    std::vector<std::vector<uint32_t>> m;
    */

    std::ifstream tableFile(inputFile);
    if(tableFile.fail()){
        std::cerr << "Could not open init file" << std::endl;
        exit(1);
    }
    load_init_table(tableFile, d, s, a, m, 2*nbPairs, offset);
    tableFile.close();

    std::vector<SobolGenerator1D<uint32_t> > sobols(m.size());
    init_sobols(sobols, d, s, a, m);

//    nbPairs = uint32_t (sobols.size()) / 2;
//    nbPairs =  4;
    std::vector<std::vector<double>> resOT(nbPairs * (octaveMax+1), std::vector<double>(nbOwens));
    std::vector<std::vector<double>> resSOT(nbPairs * (octaveMax+1), std::vector<double>(nbOwens));
    std::vector<std::vector<double>> resInte(nbPairs * (octaveMax+1), std::vector<double>(nbOwens));
    std::vector<std::vector<double>> resDisc(nbPairs * (octaveMax+1), std::vector<double>(nbOwens));

    //std::cout << "#pairIndices #nbpts #MSE #SOTCost #OTCost #Discrepancy (initFile: " << inputFile << ")" << std::endl;
    std::cout << "#pairIndices #nbpts #MSE #SOTCost (initFile: " << inputFile << ")" << std::endl;
    //std::cout << "#pairIndices #nbpts #MSE #OTCost (initFile: " << inputFile << ")" << std::endl;

    GEO::initialize();
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::set_arg("algo:predicates", "exact");
    //GEO::Process::enable_multithreading(false);
    RNG gen(seed);
    std::vector<int> seeds(nbPairs);
    for (int &v : seeds){
	    v = gen();
    }
#pragma omp parallel for
    for (int i = 0; i < nbPairs; ++i){
//       std::vector<SobolGenerator1D<uint32_t> > sobol2D(2);
//        sobol2D[0] = sobols[2 * i];
//        sobol2D[1] = sobols[2 * i + 1];
//#pragma omp parallel for
        for (uint32_t j = octaveDepart; j <= octaveMax; j += 2){
            uint32_t nbpts = 1u << j;
	    int tmpseed = seeds[i];
            resSOT[i * octaveMax + j] = testWN<2>(computeSOTCost, nbpts, nbOwens, tmpseed);
            resInte[i * octaveMax + j] = testWN<2>(computeIntegration, nbpts, nbOwens, tmpseed);
        }

#pragma omp critical (output)
        for (uint32_t j = octaveDepart, nbpts=nbptsDepart; j <= octaveMax; j += 2, nbpts *= 4) {
            for (uint32_t k = 0; k < nbOwens; ++k) {
                std::cout << i + offset/2 << ' '
                          << nbpts << ' '
                          << resInte[i * octaveMax + j][k] << ' '
                          << resSOT[i * octaveMax + j][k] << ' '
                          //<< resOT[i * octaveMax + j][k]
                          << std::endl;
            }
        }

    }

    return 0;
}
