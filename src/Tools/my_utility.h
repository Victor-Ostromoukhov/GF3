//
// Created by lpaulin on 06/02/19.
//

#ifndef SLICEDOPTIM_MY_UTILITY_H
#define SLICEDOPTIM_MY_UTILITY_H

#include <random>
#include <functional>
#include <string>

double inverseFunction(std::function<double(double)> &f, std::function<double(double)> &df, double v);

#endif //SLICEDOPTIM_MY_UTILITY_H
