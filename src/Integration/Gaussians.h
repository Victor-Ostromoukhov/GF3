//
// Created by lois on 18/09/2020.
//

#ifndef SOBOL_GAUSSIANS_H
#define SOBOL_GAUSSIANS_H

#include <cmath>
#include <vector>
#include <iostream>
#include "../Points/VecX.h"
#ifdef _MSC_VER 
#include "Tools/drand48.h" 
#endif

struct t_GaussianStruct1D {double integral; double mu[1] ; double mxCInv[1 * 1] ; };
struct t_GaussianStruct2D {double integral; double mu[2] ; double mxCInv[2 * 2] ; };
struct t_GaussianStruct3D {double integral; double mu[3] ; double mxCInv[3 * 3] ; };
struct t_GaussianStruct4D {double integral; double mu[4] ; double mxCInv[4 * 4] ; };
struct t_GaussianStruct5D {double integral; double mu[5] ; double mxCInv[5 * 5] ; };
struct t_GaussianStruct6D {double integral; double mu[6] ; double mxCInv[6 * 6] ; };
struct t_GaussianStruct8D {double integral; double mu[8] ; double mxCInv[8 * 8] ; };
struct t_GaussianStruct10D {double integral; double mu[10] ; double mxCInv[10 * 10] ; };
struct t_GaussianStruct12D {double integral; double mu[12] ; double mxCInv[12 * 12] ; };

//~ #include "../includes/Gauss2D.hpp"
extern t_GaussianStruct1D tab_Gauss1D[16384] ;
extern t_GaussianStruct2D tab_Gauss2D[16384] ;

inline double getMultivariateGaussian(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++)accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return exp(-.5 * accumulator);
}	// getMultivariateGaussian


//getMultivariateGaussian2DIsotropicCircular[{x_,y_},{mu_,mxCInv_},mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0,  Quiet[1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ]

inline double getMultivariateGaussian2DIsotropicCircular(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    double x = (pt_ND[0]-.5);
    double y = (pt_ND[1]-.5);
    if(x*x + y*y > .25) return 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++)accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return exp(-.5 * accumulator);
}	// getMultivariateGaussian2DIsotropicCircular

//getMultivariateGaussian2DIsotropicDome[{x_,y_},{mu_,mxCInv_},mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0, Quiet[Sqrt[.25 - (x-.5)^2 - (y-.5)^2] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ]

inline double getMultivariateGaussian2DIsotropicDome(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    double x = (pt_ND[0]-.5);
    double y = (pt_ND[1]-.5);
    if(x*x + y*y > .25) return 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++)accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return sqrt(.25 - x*x - y*y) * exp(-.5 * accumulator);
}	// getMultivariateGaussian2DIsotropicDome

template <int N>

inline double calculate_mse_Gaussian(const std::vector<VecX<N>>& points, const int Gaussian_mag = 5, const int nintegrands = 1024, int seed=-1) {
    double mse_accumulator = 0.;

    RNG rng;
    if (seed != -1){
        rng.seed(seed);
    } else {
        rng.seed((1u<<31u) * drand48());
    }

    int k_nintegrands = 16384 / nintegrands;
    // integration over 1K integrands over 16K
#pragma omp parallel for reduction(+:mse_accumulator)
    for (int iintegrands = 0; iintegrands < nintegrands; iintegrands ++) {
        double integral, *mu, *mxCInv;
        int integrand_index = iintegrands*k_nintegrands + floor(k_nintegrands*rng.sample_double());
        switch (N) {
            case 1 : integral = tab_Gauss1D[integrand_index].integral;
                mu = tab_Gauss1D[integrand_index].mu;
                mxCInv = tab_Gauss1D[integrand_index].mxCInv;
                break;
            case 2 : integral = tab_Gauss2D[integrand_index].integral;
                mu = tab_Gauss2D[integrand_index].mu;
                mxCInv = tab_Gauss2D[integrand_index].mxCInv;
                break;
            default: std::cout << N << "-D not implemented" << std::endl;
                exit(1);
                break;
        }
        double accumulator = 0.;
        for (const VecX<N> & point : points) {
        	accumulator += getMultivariateGaussian(N, point.data(), mu, mxCInv);
        }
        accumulator /= (double)points.size();
        double mse = (integral-accumulator)*(integral-accumulator);
        mse_accumulator += mse;
    }
    return mse_accumulator / (double)(nintegrands);
} // calculate_mse_Gaussian


#endif //SOBOL_GAUSSIANS_H
