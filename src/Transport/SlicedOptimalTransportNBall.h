#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by lpaulin on 22/03/19.
//

#ifndef SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
#define SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <omp.h>
#include "NBallRadonManager.h"
#include "../Points/VecX.h"


inline void getInverseRadonNBall(int D, int nbSamples, std::vector<double>& pos){

    pos.resize(nbSamples);
    NBallRadonManager nbrm(D);

#pragma omp parallel for shared(pos)
    for (int i = 0; i < nbSamples; ++i){
        double p = (2. * i + 1.) / (2. * nbSamples);
        pos[i] = nbrm.inverseCDF(p);
    }
}

/**
 * Compute the 1D transport cost of the from given ND \param points
 * to the distribution whose optimal 1D positions along dir are given in \param pos
 *
 * @tparam VECTYPE any vector type with * being the dot product and [] granting access to the coordinates value
 * @param points ND points
 * @param dir Projection direction
 * @param pos Sorted optimal 1D positions for target distribution along dir
 * @param pointsProject Pre allocated buffer used to sort the points's projection
 */
template <class VECTYPE>
inline double computeSlicedOTCost(const std::vector<VECTYPE>& points,
                                  const VECTYPE& dir,
                                  const std::vector<double>& pos,
                                  std::vector<double>& pointsProject
                                  ){
    double cost = 0.;

    for (size_t i = 0; i < pointsProject.size(); ++i){
        pointsProject[i] = points[i] * dir;
    }
    sort(pointsProject.begin(), pointsProject.end());

    //Place them at optimal places
    for (size_t i = 0; i < pointsProject.size(); ++i) {
        cost += (pointsProject[i] - pos[i]) * (pointsProject[i] - pos[i]);
    }
    return cost / points.size();
}

/**
 * Compute the semi discrete sliced optimal transport energy for given \param points using \param nbSlicesPerThread * \param nbThreads slices
 *
 * @tparam VECTYPE any vector type with * being the dot product and [] granting access to the coordinates value and existing method .normalize()
 * @tparam RNG any random number generator that can be used by a stl random distribution
 * @param points ND points
 * @param nbSlicesPerThread Number of slices per thread
 * @param nbThreads Number of thread to use
 * @param gen random number generator
 * @return Sliced Optimal Transport Energy
 */
template <class VECTYPE, class RNG>
inline double computeSlicedOTCost(const std::vector<VECTYPE> &points, int nbSlicesPerThread, int nbThreads, RNG &gen){

    double cost = 0.;
    int dim = points.front().dim();

    //Compute optimal 1D position for given number of points;
    int nbSamples = int(points.size());
    std::vector<double> pos(points.size());
    getInverseRadonNBall(dim, nbSamples, pos);

    //Init RNG for each thread
    std::uniform_int_distribution<int> unif;
    std::vector<std::mt19937> generators(nbThreads);
    for (int thread = 0; thread < nbThreads; ++thread){
        generators[thread].seed(unif(gen));
    }

    //Pre allocate space for projection value sorting
    std::vector<std::vector<double>> pointsProject(nbThreads, std::vector<double>(nbSamples));

    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(nbThreads); // Use 4 threads for all consecutive parallel regions
#pragma omp parallel for shared(points, pointsProject, generators) reduction(+:cost)
    for (int thread = 0; thread < nbThreads; ++thread){
        for (int slice = 0; slice < nbSlicesPerThread; ++slice){
            VECTYPE dir = randomUnitVector<VECTYPE, std::mt19937>(dim, generators[thread]);
            cost += computeSlicedOTCost(points, dir, pos, pointsProject[thread]);
        }
    }
    return cost / (nbSlicesPerThread * nbThreads);

}

/**
 * Compute the 1D discrete transport cost of the from given ND \param points
 * to the distribution whose optimal 1D positions along dir are given in \param pos
 *
 * @tparam VECTYPE any vector type with * being the dot product and [] granting access to the coordinates value
 * @param points ND points
 * @param dir Projection direction
 * @param pos Sorted optimal 1D positions for target distribution along dir
 * @param pointsProject Pre allocated buffer used to sort the points's projection
 */
template <class VECTYPE>
inline double computeDiscreteSlicedOTCost(const std::vector<VECTYPE>& points1,
                                          const std::vector<VECTYPE>& points2,
                                          const VECTYPE& dir,
                                          std::vector<double>& pointsProject1,
                                          std::vector<double>& pointsProject2
){
    double cost = 0.;

    for (size_t i = 0; i < pointsProject1.size(); ++i){
        pointsProject1[i] = points1[i] * dir;
        pointsProject2[i] = points2[i] * dir;
    }
    sort(pointsProject1.begin(), pointsProject1.end());
    sort(pointsProject2.begin(), pointsProject2.end());

    //Place them at optimal places
    for (size_t i = 0; i < pointsProject1.size(); ++i) {
        cost += (pointsProject1[i] - pointsProject2[i]) * (pointsProject1[i] - pointsProject2[i]);
    }
    return cost;
}

/**
 * Compute the discrete sliced optimal transport energy for given \param points using \param nbSlicesPerThread * \param nbThreads slices
 *
 * @tparam VECTYPE any vector type with * being the dot product and [] granting access to the coordinates value and existing method .normalize()
 * @tparam RNG any random number generator that can be used by a stl random distribution
 * @param points ND points
 * @param nbSlicesPerThread Number of slices per thread
 * @param nbThreads Number of thread to use
 * @param gen random number generator
 * @return Sliced Optimal Transport Energy
 */
template <class VECTYPE, class RNG>
inline double computeDiscreteSlicedOTCost(const std::vector<VECTYPE>& points1, const std::vector<VECTYPE>& points2, int nbSlicesPerThread, int nbThreads, RNG& gen){

    double cost = 0.;
    int dim = points1.front().dim();

    int nbSamples = int(points1.size());

    //Init RNG for each thread
    std::uniform_int_distribution<int> unif;
    std::vector<std::mt19937> generators(nbThreads);
    for (int thread = 0; thread < nbThreads; ++thread){
        generators[thread].seed(unif(gen));
    }

    //Pre allocate space for projection value sorting
    std::vector<std::vector<double>> pointsProject1(nbThreads, std::vector<double>(nbSamples));
    std::vector<std::vector<double>> pointsProject2(nbThreads, std::vector<double>(nbSamples));

#pragma omp parallel for shared(points1, points2, pointsProject1, pointsProject2, generators) reduction(+:cost)
    for (int thread = 0; thread < nbThreads; ++thread){
        for (int slice = 0; slice < nbSlicesPerThread; ++slice){
            VECTYPE dir = randomUnitVector<VECTYPE, std::mt19937>(dim, generators[thread]);
            cost += computeDiscreteSlicedOTCost(points1, points2, dir, pointsProject1[thread], pointsProject2[thread]);
        }
    }
    return cost / (nbSlicesPerThread * nbThreads);

}

#endif //SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H

#pragma clang diagnostic pop