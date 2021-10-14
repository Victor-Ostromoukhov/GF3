//
// Created by lois on 18/09/2020.
//

#ifndef SOBOL_GAUSSIANS_H
#define SOBOL_GAUSSIANS_H

#include <cmath>
#include <vector>
#include <iostream>
#include "../Points/VecX.h"

struct t_Heaviside2D {long double integral; long double muDiscotinuity[2]; long double normal[2]; };
struct t_Heaviside3D {long double integral; long double muDiscotinuity[3]; long double normal[3]; };
struct t_Heaviside4D {long double integral; long double muDiscotinuity[4]; long double normal[4]; };
struct t_Heaviside6D {long double integral; long double muDiscotinuity[6]; long double normal[6]; };

struct t_Gauss2D {long double integral; long double mu[2] ; long double mxCInv[2 * 2] ; };
struct t_Gauss3D {long double integral; long double mu[3] ; long double mxCInv[3 * 3] ; };
struct t_Gauss4D {long double integral; long double mu[4] ; long double mxCInv[4 * 4] ; };
struct t_Gauss6D {long double integral; long double mu[6] ; long double mxCInv[6 * 6] ; };

struct t_GaussianStruct2D {
	long double integral;
	long double mu[2];
	long double mxCInv[2 * 2];
};
struct t_HeaviBell2D {
	long double integral;
	long double muDiscotinuity[2];
	long double normal[2];
};
struct t_HeaviCont2D {
	long double integral;
	long double muDiscotinuity[2];
	long double normal[2];
	long double mu[2];
	long double mxCInv[2 * 2];
};
struct t_HeaviGauss2D {
	long double integral;
	long double muDiscotinuity[2];
	long double normal[2];
	long double mu[2];
	long double mxCInv[2 * 2];
};


#include "../includes/Heaviside2D.hpp"
#include "../includes/Heaviside3D.hpp"
#include "../includes/Heaviside4D.hpp"
#include "../includes/Heaviside6D.hpp"

#include "../includes/Gauss2D.hpp"
#include "../includes/Gauss3D.hpp"
#include "../includes/Gauss4D.hpp"
#include "../includes/Gauss6D.hpp"

#include "../includes/Smooth2D.hpp"
#include "../includes/SmoothBump2D.hpp"
#include "../includes/Cont2D.hpp"
#include "../includes/GaussIso2D.hpp"
#include "../includes/HeaviBell2D.hpp"
#include "../includes/HeaviGauss2D.hpp"
#include "../includes/HeaviCont2D.hpp"

long double getMultivariateGaussian(const int nDims, const double *pt_ND, const long double *mu, const long double *mxCInv) {
    long double accumulator = 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    long double res = exp(-.5 * accumulator);
    return res;
}	// getMultivariateGaussian

#define NDIMS 2
// getMultivariateSmooth2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := With[{r2 = (x-.5)^2 + (y-.5)^2 }, If[r2 > .25, 0,  54.5982 Quiet[Exp[-1/((.5)^4 - r2^2)] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ] (* bump function https://en.m.wikipedia.org/wiki/Bump_function *)]
inline long double getSmooth2D(const double *pt_ND, const long double *mu, const long double *mxCInv) {
    long double accumulator = 0.;
    long double x = (pt_ND[0]-.5);
    long double y = (pt_ND[1]-.5);
    long double r2 = x*x + y*y;
    if(r2 >= .25) return 0.;
    for (int row = 0; row < NDIMS; row++) {
        for (int col = 0; col < NDIMS; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*NDIMS+row];
    }
    return 8886110.520507872*exp(-1./(.0625 - r2*r2 ) ) * exp(-.5 * accumulator);
}	// getSmooth2D

// getMultivariateSmoothBump2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := With[{r2 = (x-.5)^2 + (y-.5)^2 }, If[r2 > .25, 0,  54.5982 Quiet[Exp[-1/((.5)^2 - r2)] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ] (* bump function https://en.m.wikipedia.org/wiki/Bump_function *)]
inline long double getSmoothBump2D(const double *pt_ND, const long double *mu, const long double *mxCInv) {
    long double accumulator = 0.;
    long double x = (pt_ND[0]-.5);
    long double y = (pt_ND[1]-.5);
    long double r2 = x*x + y*y;
    if(r2 >= .25) return 0.;
    for (int row = 0; row < NDIMS; row++) {
        for (int col = 0; col < NDIMS; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*NDIMS+row];
    }
    return 54.5982*exp(-1./(.25 - r2 ) ) * exp(-.5 * accumulator);
}	// getSmoothBump2D

//getCont2D[{x_,y_},{mu_,mxCInv_},mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0, Quiet[Sqrt[.25 - (x-.5)^2 - (y-.5)^2] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ]
inline long double getCont2D(const double *pt_ND, const long double *mu, const long double *mxCInv) {
    long double accumulator = 0.;
    long double x = (pt_ND[0]-.5);
    long double y = (pt_ND[1]-.5);
    long double r2 = x*x + y*y;
    if(r2 >= .25) return 0.;
    for (int row = 0; row < NDIMS; row++) {
        for (int col = 0; col < NDIMS; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*NDIMS+row];
    }
    return sqrt(.25 - r2) * exp(-.5 * accumulator);
}	// getCont2D

#define NINTEGRANS	1024

template <typename Vec>
inline long double calculate_mse(const std::vector<Vec>& points, const int integrandType = 1, const int nintegrands = 1024) {
    long double mse_accumulator = 0.;
    int k_nintegrands = NINTEGRANS / nintegrands;
    int nDims = points[0].dim();
    // integration over 1K integrands over 16K
//#pragma omp parallel for schedule(dynamic) reduction(+:mse_accumulator)
#pragma omp parallel for reduction(+:mse_accumulator)
    for (int iintegrands = 0; iintegrands < nintegrands; iintegrands ++) {
        long double integral, *mu, *mxCInv, *muDiscotinuity, *normal;
        int integrand_index = iintegrands*k_nintegrands + floor(k_nintegrands*drand48());
        long double accumulator = 0.;
        switch (integrandType) {
        	case 1 : // Heaviside2D
            	// getHeaviside2D[{x_,y_},{muDiscotinuity_,normVector_},mulFactor_:1] := If[({x,y}-muDiscotinuity).normVector > 0, 1, 0]
        		switch(nDims) {
            	case 2 : // Heaviside2D
                  	integral = tab_Heaviside2D[integrand_index].integral;
                  	muDiscotinuity = tab_Heaviside2D[integrand_index].muDiscotinuity;
                  	normal = tab_Heaviside2D[integrand_index].normal;
                  	break;
            	case 3 : // Heaviside3D
                  	integral = tab_Heaviside3D[integrand_index].integral;
                  	muDiscotinuity = tab_Heaviside3D[integrand_index].muDiscotinuity;
                  	normal = tab_Heaviside3D[integrand_index].normal;
            		break;
            	case 4 : // Heaviside4D
                  	integral = tab_Heaviside4D[integrand_index].integral;
                  	muDiscotinuity = tab_Heaviside4D[integrand_index].muDiscotinuity;
                  	normal = tab_Heaviside4D[integrand_index].normal;
                  	break;
            	case 6 : // Heaviside6D
                  	integral = tab_Heaviside6D[integrand_index].integral;
                  	muDiscotinuity = tab_Heaviside6D[integrand_index].muDiscotinuity;
                  	normal = tab_Heaviside6D[integrand_index].normal;
            		break;
                default: std::cerr << "nDims=" << nDims << ": not implemented." << std::endl;
                      exit(1);
                      break;
        		}
                for (int ipt = 0; ipt < points.size() ; ipt++) {
                	long double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < nDims ; idim++) {
                       	long double integrand_mu_componenet = muDiscotinuity[idim];
                       	long double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     if (dotProduct > 0 ) accumulator += 1.;
                }
                break;
            case 2 :	// Gauss (slightly anisotropic)
               	// getMultivariate2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := Quiet[ 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)] ]
        		switch(nDims) {
            	case 2 : // Heaviside2D
                	integral = tab_Gauss2D[integrand_index].integral;
                    mu = tab_Gauss2D[integrand_index].mu;
                    mxCInv = tab_Gauss2D[integrand_index].mxCInv;
                  	break;
            	case 3 : // Heaviside3D
                	integral = tab_Gauss3D[integrand_index].integral;
                    mu = tab_Gauss3D[integrand_index].mu;
                    mxCInv = tab_Gauss3D[integrand_index].mxCInv;
            		break;
            	case 4 : // Heaviside4D
                	integral = tab_Gauss4D[integrand_index].integral;
                    mu = tab_Gauss4D[integrand_index].mu;
                    mxCInv = tab_Gauss4D[integrand_index].mxCInv;
                  	break;
            	case 6 : // Heaviside6D
                	integral = tab_Gauss6D[integrand_index].integral;
                    mu = tab_Gauss6D[integrand_index].mu;
                    mxCInv = tab_Gauss6D[integrand_index].mxCInv;
            		break;
                default: std::cerr << "nDims=" << nDims << ": not implemented." << std::endl;
                      exit(1);
                      break;
        		}
                for (const Vec & point : points) {
                	accumulator += getMultivariateGaussian(nDims, point.data(), mu, mxCInv);
                }
                break;
            case 3 :	// GaussIso (isotropic)
            	// getMultivariate2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := Quiet[ 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)] ]
            	integral = tab_GaussIso[integrand_index].integral;
                mu = tab_GaussIso[integrand_index].mu;
                mxCInv = tab_GaussIso[integrand_index].mxCInv;
                for (const Vec & point : points) {
                	accumulator += getMultivariateGaussian(2, point.data(), mu, mxCInv);
                }
                break;
        	case 4 :	// Cont2D
               	// getMultivariateContND[{x_,y_},{mu_,mxCInv_},mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0, Quiet[Sqrt[.25 - (x-.5)^2 - (y-.5)^2] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ]
            	integral = tab_Cont2D[integrand_index].integral;
                mu = tab_Cont2D[integrand_index].mu;
                mxCInv = tab_Cont2D[integrand_index].mxCInv;
                for (const Vec & point : points) {
                	accumulator += getCont2D(point.data(), mu, mxCInv);
                }
            break;
            case 5 : // HeaviBell2D
              	// getHeaviBell2D[{x_,y_},muDiscotinuity_,normVector_,mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0, If[({x,y}-muDiscotinuity).normVector > 0, Sqrt[.25 - (x-.5)^2 - (y-.5)^2], 0] ]
				integral = tab_HeaviBell2D[integrand_index].integral;
				muDiscotinuity = tab_HeaviBell2D[integrand_index].muDiscotinuity;
				normal = tab_HeaviBell2D[integrand_index].normal;
				for (int ipt = 0; ipt < points.size() ; ipt++) {
					long double dotProduct = 0;
					Vec pt = points[ipt];
	                   for (int idim = 0; idim < pt.dim() ; idim++) {
	                     	long double integrand_mu_componenet = muDiscotinuity[idim];
	                     	long double integrand_normal_componenet = normal[idim];
	                     	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
	                   }
                   long double xhalf = pt[0] - .5;
                   long double yhalf = pt[1] - .5;
                   long double xhalf2 = xhalf*xhalf;
                   long double yhalf2 = yhalf*yhalf;
                   if (dotProduct > 0 && xhalf2 + yhalf2 < .25) accumulator += sqrt(.25 - xhalf2 - yhalf2);
	              }
	              break;
            case 6 : // HeaviCont
               	// getHeaviCont2D[{x_,y_},muDiscotinuity_,normVector_,mu_,mxCInv_,mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0, If[({x,y}-muDiscotinuity).normVector > 0, Sqrt[.25 - (x-.5)^2 - (y-.5)^2] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)], 0] ]
            	integral = tab_HeaviCont2D[integrand_index].integral;
            	muDiscotinuity = tab_HeaviCont2D[integrand_index].muDiscotinuity;
            	normal = tab_HeaviCont2D[integrand_index].normal;
                mu = tab_HeaviCont2D[integrand_index].mu;
                mxCInv = tab_HeaviCont2D[integrand_index].mxCInv;
                for (int ipt = 0; ipt < points.size() ; ipt++) {
                	long double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < pt.dim() ; idim++) {
                       	long double integrand_mu_componenet = muDiscotinuity[idim];
                       	long double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     long double xhalf = pt[0] - .5;
                     long double yhalf = pt[1] - .5;
                     long double xhalf2 = xhalf*xhalf;
                     long double yhalf2 = yhalf*yhalf;
                     if (dotProduct > 0 && xhalf2 + yhalf2 < .25) accumulator += getCont2D(pt.data(), mu, mxCInv);
                }
                break;
            case 7 : // HeaviGauss2D
               	// getHeaviGauss2D[{x_,y_},{muDiscotinuity_,normVector_,mu_,mxCInv_},mulFactor_:1] := If[({x,y}-muDiscotinuity).normVector > 0, 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)], 0]
            	integral = tab_HeaviGauss2D[integrand_index].integral;
            	muDiscotinuity = tab_HeaviGauss2D[integrand_index].muDiscotinuity;
            	normal = tab_HeaviGauss2D[integrand_index].normal;
                mu = tab_HeaviGauss2D[integrand_index].mu;
                mxCInv = tab_HeaviGauss2D[integrand_index].mxCInv;
                for (int ipt = 0; ipt < points.size() ; ipt++) {
                	long double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < pt.dim() ; idim++) {
                       	long double integrand_mu_componenet = muDiscotinuity[idim];
                       	long double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     if (dotProduct > 0) accumulator += getMultivariateGaussian(2, points[ipt].data(), mu, mxCInv);
                }
                break;
            case 8 :	// SmoothBump2D using bump function exp(-1/(1-r^4))  https://en.m.wikipedia.org/wiki/Bump_function
            	// getMultivariateSmoothBump2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := With[{r2 = (x-.5)^2 + (y-.5)^2 }, If[r2 > .25, 0,  54.5982 Quiet[Exp[-1/((.5)^2 - r2)] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ] (* bump function https://en.m.wikipedia.org/wiki/Bump_function *)]
            	integral = tab_SmoothBump2D[integrand_index].integral;
                mu = tab_SmoothBump2D[integrand_index].mu;
                mxCInv = tab_SmoothBump2D[integrand_index].mxCInv;
                for (const Vec   & point : points) {
                	accumulator += getSmoothBump2D(point.data(), mu, mxCInv);
                }
                break;
            case 9 :	// Smooth2D using bump function exp(-1/(1-r^4))  https://en.m.wikipedia.org/wiki/Bump_function
            	// getMultivariateSmooth2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := With[{r2 = (x-.5)^2 + (y-.5)^2 }, If[r2 > .25, 0,  54.5982 Quiet[Exp[-1/((.5)^4 - r2^2)] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ] (* bump function https://en.m.wikipedia.org/wiki/Bump_function *)]
            	integral = tab_Smooth2D[integrand_index].integral;
                mu = tab_Smooth2D[integrand_index].mu;
                mxCInv = tab_Smooth2D[integrand_index].mxCInv;
                for (const Vec & point : points) {
                	accumulator += getSmooth2D(point.data(), mu, mxCInv);
                }
                break;
          default: std::cout << "integrandType " << integrandType << " not implemented." << std::endl;
                exit(1);
                break;
        }
        accumulator /= (long double)points.size();
        long double mse = (integral-accumulator)*(integral-accumulator);
        mse_accumulator += mse;
    }
    return mse_accumulator / (long double)(nintegrands);
} // calculate_mse


#endif //SOBOL_GAUSSIANS_H
