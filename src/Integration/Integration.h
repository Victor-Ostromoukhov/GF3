//
// Created by lois on 18/09/2020.
//

#ifndef SOBOL_GAUSSIANS_H
#define SOBOL_GAUSSIANS_H

#include <cmath>
#include <vector>
#include <iostream>
#include "../Points/VecX.h"
#include "../Samplers/Random.h"
#ifdef _MSC_VER
#include "Tools/drand48.h" 
#endif

struct t_Heaviside2D {double integral; double muDiscotinuity[2]; double normal[2]; };
struct t_Heaviside3D {double integral; double muDiscotinuity[3]; double normal[3]; };
struct t_Heaviside4D {double integral; double muDiscotinuity[4]; double normal[4]; };
struct t_Heaviside5D {double integral; double muDiscotinuity[5]; double normal[5]; };
struct t_Heaviside6D {double integral; double muDiscotinuity[6]; double normal[6]; };

struct t_Gauss1D {double integral; double mu[1] ; double mxCInv[1 * 1] ; };
struct t_Gauss2D {double integral; double mu[2] ; double mxCInv[2 * 2] ; };
struct t_Gauss3D {double integral; double mu[3] ; double mxCInv[3 * 3] ; };
struct t_Gauss4D {double integral; double mu[4] ; double mxCInv[4 * 4] ; };
struct t_Gauss5D {double integral; double mu[5] ; double mxCInv[5 * 5] ; };
struct t_Gauss6D {double integral; double mu[6] ; double mxCInv[6 * 6] ; };
struct t_Gauss8D {double integral; double mu[8] ; double mxCInv[8 * 8] ; };
struct t_Gauss10D {double integral; double mu[10] ; double mxCInv[10 * 10] ; };
struct t_Gauss12D {double integral; double mu[12] ; double mxCInv[12 * 12] ; };

struct t_GaussianStruct1D {
    double integral;
    double mu[1];
    double mxCInv[1 * 1];
};
struct t_GaussianStruct2D {
    double integral;
    double mu[2];
    double mxCInv[2 * 2];
};
struct t_HeaviBell2D {
	double integral;
	double muDiscotinuity[2];
	double normal[2];
};
struct t_HeaviCont2D {
	double integral;
	double muDiscotinuity[2];
	double normal[2];
	double mu[2];
	double mxCInv[2 * 2];
};
struct t_HeaviGauss2D {
	double integral;
	double muDiscotinuity[2];
	double normal[2];
	double mu[2];
	double mxCInv[2 * 2];
};


//~ #include "../includes/Heaviside2D.hpp"
extern t_Heaviside2D tab_Heaviside2D[16384] ;
//~ #include "../includes/Heaviside3D.hpp"
extern t_Heaviside3D tab_Heaviside3D[16384] ;
//~ #include "../includes/Heaviside4D.hpp"
extern t_Heaviside4D tab_Heaviside4D[16384] ;
//~ #include "../includes/Heaviside5D.hpp"
extern t_Heaviside5D tab_Heaviside5D[16384] ;
//~ #include "../includes/Heaviside6D.hpp"
extern t_Heaviside6D tab_Heaviside6D[16384] ;

//~ #include "../includes/Gauss2D.hpp"
extern t_Gauss1D tab_Gauss1D[16384] ;
extern t_GaussianStruct2D tab_Gauss2D[16384] ;
//~ #include "../includes/Gauss3D.hpp"
extern t_Gauss3D tab_Gauss3D[16384] ;
//~ #include "../includes/Gauss4D.hpp"
extern t_Gauss4D tab_Gauss4D[16384] ;
//~ #include "../includes/Gauss5D.hpp"
extern t_Gauss5D tab_Gauss5D[16384] ;
//~ #include "../includes/Gauss6D.hpp"
extern t_Gauss6D tab_Gauss6D[16384] ;
extern t_Gauss8D tab_Gauss8D[16384] ;
extern t_Gauss10D tab_Gauss10D[16384] ;
extern t_Gauss12D tab_Gauss12D[16384] ;

//~ #include "../includes/Smooth2D.hpp"
extern t_Gauss2D tab_Smooth2D[16384] ;
//~ #include "../includes/Smooth3D.hpp"
extern t_Gauss3D tab_Smooth3D[16384] ;
//~ #include "../includes/Smooth4D.hpp"
extern t_Gauss4D tab_Smooth4D[16384] ;
//~ #include "../includes/Cont2D.hpp"
extern t_GaussianStruct2D tab_Cont2D[16384] ;
//~ #include "../includes/GaussIso2D.hpp"
extern t_GaussianStruct2D tab_GaussIso[16384] ;
//~ #include "../includes/HeaviBell2D.hpp"
extern t_HeaviBell2D tab_HeaviBell2D[16384] ;
//~ #include "../includes/HeaviGauss2D.hpp"
extern t_HeaviGauss2D tab_HeaviGauss2D[16384] ;
//~ #include "../includes/HeaviCont2D.hpp"
extern t_HeaviCont2D tab_HeaviCont2D[16384] ;

inline double getMultivariateGaussian(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return exp(-.5 * accumulator);
}	// getMultivariateGaussian

inline double getSmoothND(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0., r2 = 0.;
    for (int idim = 0; idim < nDims; idim++)
    	r2 += ((pt_ND[idim]-.5)*(pt_ND[idim]-.5));
    if(r2 >= .25) return 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return 54.5982*exp(-1./(.25 - r2 ) ) * exp(-.5 * accumulator);
}	// getSmoothND

inline double getSmoothPow4ND(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0., r2 = 0.;
    for (int idim = 0; idim < nDims; idim++)
    	r2 += ((pt_ND[idim]-.5)*(pt_ND[idim]-.5));
    if(r2 >= .25) return 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return 8886110.520507872*exp(-1./(.0625 - r2*r2 ) ) * exp(-.5 * accumulator);
}	// getSmoothPow4ND

#define NDIMS 2
//getCont2D[{x_,y_},{mu_,mxCInv_},mulFactor_:1] := If[(x-.5)^2 + (y-.5)^2 > .25, 0, Quiet[Sqrt[.25 - (x-.5)^2 - (y-.5)^2] 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)]] ]
inline double getCont2D(const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    double x = (pt_ND[0]-.5);
    double y = (pt_ND[1]-.5);
    double r2 = x*x + y*y;
    if(r2 >= .25) return 0.;
    for (int row = 0; row < NDIMS; row++) {
        for (int col = 0; col < NDIMS; col++) accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*NDIMS+row];
    }
    return sqrt(.25 - r2) * exp(-.5 * accumulator);
}	// getCont2D

template <typename Vec>
inline double calculate_mse(const std::vector<Vec>& points,
                            const int integrandType = 1,
                            const int nintegrands = 1024,
                            const int seed=-1)
{
    RNG rng;
    if (seed != -1){
        rng.seed(seed);
    } else {
        rng.seed((1u<<31u) * drand48());
    }

    double mse_accumulator = 0.;
    int k_nintegrands = 16384 / nintegrands;
    int nDims = points[0].dim();
    // integration over 1K integrands over 16K
//#pragma omp parallel for schedule(dynamic) reduction(+:mse_accumulator)
#pragma omp parallel for reduction(+:mse_accumulator)
    for (int iintegrands = 0; iintegrands < nintegrands; iintegrands ++) {
        double integral, *mu, *mxCInv, *muDiscotinuity, *normal;
        int integrand_index;
        integrand_index = iintegrands*k_nintegrands + floor(k_nintegrands * rng.sample_double());// drand48());
        if(nintegrands == 1)
        	integrand_index = 1;
        else
            integrand_index = iintegrands*k_nintegrands + floor(k_nintegrands * rng.sample_double());// drand48());
        double accumulator = 0.;
        switch (integrandType) {
        	case 1 : // HeavisideND
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
            	case 5 : // Heaviside5D
                  	integral = tab_Heaviside5D[integrand_index].integral;
                  	muDiscotinuity = tab_Heaviside5D[integrand_index].muDiscotinuity;
                  	normal = tab_Heaviside5D[integrand_index].normal;
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
                	double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < nDims ; idim++) {
                       	double integrand_mu_componenet = muDiscotinuity[idim];
                       	double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     if (dotProduct > 0 ) accumulator += 1.;
                }
                break;
            case 2 :	// Gauss (slightly anisotropic)
               	// getMultivariate2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := Quiet[ 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)] ]
        		switch(nDims) {
            	case 1 : // Gauss 1D
                	integral = tab_Gauss1D[integrand_index].integral;
                    mu = tab_Gauss1D[integrand_index].mu;
                    mxCInv = tab_Gauss1D[integrand_index].mxCInv;
                  	break;
            	case 2 : // Gauss 2D
                	integral = tab_Gauss2D[integrand_index].integral;
                    mu = tab_Gauss2D[integrand_index].mu;
                    mxCInv = tab_Gauss2D[integrand_index].mxCInv;
                  	break;
            	case 3 : // Gauss 3D
                	integral = tab_Gauss3D[integrand_index].integral;
                    mu = tab_Gauss3D[integrand_index].mu;
                    mxCInv = tab_Gauss3D[integrand_index].mxCInv;
            		break;
            	case 4 : // Gauss 4D
                	integral = tab_Gauss4D[integrand_index].integral;
                    mu = tab_Gauss4D[integrand_index].mu;
                    mxCInv = tab_Gauss4D[integrand_index].mxCInv;
                  	break;
            	case 5 : // Gauss 5D
                	integral = tab_Gauss5D[integrand_index].integral;
                    mu = tab_Gauss5D[integrand_index].mu;
                    mxCInv = tab_Gauss5D[integrand_index].mxCInv;
            		break;
            	case 6 : // Gauss 6D
                	integral = tab_Gauss6D[integrand_index].integral;
                    mu = tab_Gauss6D[integrand_index].mu;
                    mxCInv = tab_Gauss6D[integrand_index].mxCInv;
            		break;
            	case 8 : // Gauss 8D
                	integral = tab_Gauss8D[integrand_index].integral;
                    mu = tab_Gauss8D[integrand_index].mu;
                    mxCInv = tab_Gauss8D[integrand_index].mxCInv;
            		break;
            	case 10 : // Gauss 10D
                	integral = tab_Gauss10D[integrand_index].integral;
                    mu = tab_Gauss10D[integrand_index].mu;
                    mxCInv = tab_Gauss10D[integrand_index].mxCInv;
            		break;
            	case 12 : // Gauss 12D
                	integral = tab_Gauss12D[integrand_index].integral;
                    mu = tab_Gauss12D[integrand_index].mu;
                    mxCInv = tab_Gauss12D[integrand_index].mxCInv;
            		break;
                default: std::cerr << "nDims=" << nDims << ": not implemented." << std::endl;
                      exit(1);
                      break;
        		}
                for (const Vec & point : points) {
                	accumulator += getMultivariateGaussian(nDims, point.data(), mu, mxCInv);
                }
                break;
        	case 3 :	// Smooth2D == SmoothBump2D
                	// getMultivariateSmoothND[pt_, {mu_,mxCInv_}, mulFactor_:1] := With[{r2 = euclidlen2[pt-.5]},If[r2 >= .25, 0,  8886110.520507872 Quiet[Exp[-1/((.5)^4 - r2^2)] 1./mulFactor getMultivariateND[pt,{mu,mxCInv}] ] ] ]
            		switch(nDims) {
                	case 2 : // Smooth 2D
                    	integral = tab_Smooth2D[integrand_index].integral;
                        mu = tab_Smooth2D[integrand_index].mu;
                        mxCInv = tab_Smooth2D[integrand_index].mxCInv;
                      	break;
                	case 3 : // Smooth 3D
                    	integral = tab_Smooth3D[integrand_index].integral;
                        mu = tab_Smooth3D[integrand_index].mu;
                        mxCInv = tab_Smooth3D[integrand_index].mxCInv;
                		break;
                	case 4 : // Smooth 4D
                    	integral = tab_Smooth4D[integrand_index].integral;
                        mu = tab_Smooth4D[integrand_index].mu;
                        mxCInv = tab_Smooth4D[integrand_index].mxCInv;
                      	break;
                    default: std::cerr << "nDims=" << nDims << ": not implemented." << std::endl;
                          exit(1);
                          break;
            		}
                    for (const Vec & point : points) {
                    	accumulator += getSmoothPow4ND(nDims, point.data(), mu, mxCInv);
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
					double dotProduct = 0;
					Vec pt = points[ipt];
	                   for (int idim = 0; idim < pt.dim() ; idim++) {
	                     	double integrand_mu_componenet = muDiscotinuity[idim];
	                     	double integrand_normal_componenet = normal[idim];
	                     	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
	                   }
                   double xhalf = pt[0] - .5;
                   double yhalf = pt[1] - .5;
                   double xhalf2 = xhalf*xhalf;
                   double yhalf2 = yhalf*yhalf;
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
                	double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < pt.dim() ; idim++) {
                       	double integrand_mu_componenet = muDiscotinuity[idim];
                       	double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     double xhalf = pt[0] - .5;
                     double yhalf = pt[1] - .5;
                     double xhalf2 = xhalf*xhalf;
                     double yhalf2 = yhalf*yhalf;
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
                	double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < pt.dim() ; idim++) {
                       	double integrand_mu_componenet = muDiscotinuity[idim];
                       	double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     if (dotProduct > 0) accumulator += getMultivariateGaussian(2, points[ipt].data(), mu, mxCInv);
                }
                break;
          default: std::cout << "integrandType " << integrandType << " not implemented." << std::endl;
                exit(1);
                break;
        }
        accumulator /= (double)points.size();
        double mse = (integral-accumulator)*(integral-accumulator);
        mse_accumulator += mse;
    }
    return mse_accumulator / (double)(nintegrands);
} // calculate_mse


#endif //SOBOL_GAUSSIANS_H
