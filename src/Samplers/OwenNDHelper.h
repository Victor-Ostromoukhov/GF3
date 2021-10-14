//
// Created by lois on 20/11/2020.
//

#ifndef SOBOL_OWENNDHELPER_H
#define SOBOL_OWENNDHELPER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "SobolGenerator1D.h"
#include "OwenScrambling.h"

std::vector<uint32_t> loadMapping(const std::string& filename){
    std::vector<std::pair<uint32_t, uint32_t> > affectations;
    std::ifstream in(filename);
    uint32_t a,b;
    uint32_t indMax = 0;
    char c;
    std::string dump;
    while (in >> c){
        if (c == '#'){
            std::getline(in, dump);
            continue;
        }
        in.putback(c);
        in >> a >> b;
        indMax = std::max(std::max(a, b), indMax);
        affectations.emplace_back(a,b);
    }
    std::vector<uint32_t>mapping(indMax + 1);
    for(int i = 0; i < mapping.size(); ++i){
        mapping[i] = i;
    }
    for(std::pair<uint32_t, uint32_t> &affectation : affectations){
        mapping[affectation.first] = affectation.second;
    }

    return mapping;
}

template<typename integer>
class OwenNDHelper {
public:
    std::vector<SobolGenerator1D<integer>> sobols;
    std::vector<uint32_t > seeds;
    std::vector<uint32_t > mapping;

    OwenNDHelper(std::vector<SobolGenerator1D<integer>> sobols,
                 std::vector<uint32_t > seeds,
                 std::vector<uint32_t > mapping) : sobols(sobols), seeds(seeds), mapping(mapping){}

    OwenNDHelper() = delete;

    double operator()(int dim, int idPoint, int owen_tree_depth = 8 * sizeof(integer)){
        const SobolGenerator1D<integer> &sobol = sobols[mapping[dim]];
        return OwenScrambling(sobol.getSobolInt(idPoint), seeds[mapping[dim]], owen_tree_depth);
    }

};


#endif //SOBOL_OWENNDHELPER_H
