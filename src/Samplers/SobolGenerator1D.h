// SobolGenerator1D.h
// Created by lpaulin on 05/05/20.
//

#ifndef MYSOBOL_SOBOLGENERATOR1D_H
#define MYSOBOL_SOBOLGENERATOR1D_H

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <bitset>
#include <climits>
#include <limits>
#include <sstream>
#include <fstream>

#include "OwenScrambling.h"

#ifdef _MSC_VER 
typedef unsigned int uint;
#endif

template <typename uint32_t>
inline uint32_t ReverseBits(uint32_t n);

template <>
inline uint32_t ReverseBits<uint32_t >(uint32_t n) {
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    return n;
}

template <>
inline uint64_t ReverseBits<uint64_t >(uint64_t n) {
    n = (n << 32) | (n >> 32);
    n = ((n & 0x0000ffff0000ffffULL) << 16) | ((n & 0xffff0000ffff0000ULL) >> 16);
    n = ((n & 0x00ff00ff00ff00ffULL) << 8) | ((n & 0xff00ff00ff00ff00ULL) >> 8);
    n = ((n & 0x0f0f0f0f0f0f0f0fULL) << 4) | ((n & 0xf0f0f0f0f0f0f0f0ULL) >> 4);
    n = ((n & 0x3333333333333333ULL) << 2) | ((n & 0xccccccccccccccccULL) >> 2);
    n = ((n & 0x5555555555555555ULL) << 1) | ((n & 0xaaaaaaaaaaaaaaaaULL) >> 1);
    return n;
}

template <typename uint32_t>
class SobolGenerator1D{

public:
    //dimension of Sobol polynomial
    uint32_t s;
    uint32_t d;

    //binary Sobol polynomial
    uint32_t a;

    //Sobol matrice
    std::array<uint32_t, CHAR_BIT*sizeof(uint32_t)> m;

    uint32_t owen_seeds_per_dim;

    inline SobolGenerator1D(){}

    inline SobolGenerator1D(uint32_t _d, uint32_t _s, uint32_t _a, const std::vector<uint32_t>& _m){
        init1D(_d, _s, _a, _m);
    }

    inline void generateMatrix(){
        //fill each matrix rank
//    	std::cout << "generateMatrix for dim " << d << std::endl;
    	if (s == 0) {	// van der Corput
    		for (uint32_t k = s; k < m.size(); ++k)
    			m[k] = 1;
    	} else {
            for (uint32_t k = s; k < m.size(); ++k) {
                m[k] = 0;
                for (int i = 1; i < s; ++i) {
                    // `akj` stores aj_k, note that the polynomial rep is reversed
                    // `pw2` stores 2^k
                    const uint32_t ai = (a >> (s - 1 - i)) & 1;
                    const uint32_t pw2 = (1 << i);
                    m[k] ^= ai * pw2 * m[k - i];
                }
                m[k] ^= (1 << s) * m[k - s];
                m[k] ^= m[k - s];
            }
    	}
    }


    inline void init1D(uint32_t _d, uint32_t _s, uint32_t _a, const std::vector<uint32_t>& _m, const uint32_t seed = 13374269, const bool dbg_flag = false) {
        d = _d;
        s = _s;
        a = _a;
        for (int i = 0; i < _m.size(); ++i){
            m[i] = _m[i];
        }
        generateMatrix();
//        if (dbg_flag && d < 100) {
//        	std::cout << "init1D dsa: " << d << " " << s << " " << a << "\t### ";
//        	for(int j = 0; j < m.size(); j++) {
//        		std::cout << m[j] << " ";
//        	}
//        	std::cout << std::endl;
//        }
    }	// init1D

    inline uint32_t getSobolInt(uint32_t n) const{
//        if(d == 0){	// van der Corput -> now treated by proper init in generateMatrix()
//            return ReverseBits<uint32_t>(n);
//        }
        uint32_t res = 0;
        for (int i = 0; i < CHAR_BIT * sizeof(uint32_t); ++i){
            res <<= 1;
            res ^= (n & 1) * m[i];
            n /= 2;
        }
        return res;
    }	// getSobolInt

    inline double getSobolDouble(uint32_t n) const{

        uint32_t res = getSobolInt(n);
        return double(res) / std::numeric_limits<uint32_t>::max();
    }

    inline friend std::ostream& operator<<(std::ostream& out, const SobolGenerator1D<uint32_t>& sobol){

        /*out << "a: ";
        uint32_t cpa = sobol.a;
        for (int i = 0; i < CHAR_BIT * sizeof(uint32_t); ++i){
            out << (cpa & 1);
            cpa >>= 1;
        }
        out << std::endl;
        out << "m:" << std::endl;

        for(int i = 0; i < sobol.m.size(); ++i){
            out << std::bitset<sizeof(uint32_t) * CHAR_BIT>(sobol.m[i] << (sizeof(uint32_t) * CHAR_BIT - 1 - i)) << std::endl;
        }*/

        out << sobol.d << "\t" << sobol.s << "\t" << sobol.a  << "\t";
        for (uint i = 0; i < sobol.s; ++i) {
            out << sobol.m[i] << " ";
        }

        return out;
    }

};


inline bool readOneInitVector(std::istream& in, uint32_t& d, uint32_t& s, uint32_t& a, std::vector<uint32_t>& m){

    std::string line;
    m.clear();
    //Read a line and skip if it's a header
    bool test;
    do {
        test = bool(getline(in, line));
    } while (test && line[0] == 'd');

    if(test){
        std::istringstream lineStream(line);

        lineStream >> d  >> s >> a;

        uint32_t mi;
        while(lineStream >> mi){
            m.push_back(mi);
        }
        return true;
    } else {
        return false;
    }
}

inline void load_init_table(std::istream& in,
                            std::vector<uint32_t>& d,
                            std::vector<uint32_t>& s,
                            std::vector<uint32_t>& a,
                            std::vector<std::vector<uint32_t>>& m,
                            uint32_t nbMax = -1,
                            uint32_t offset = 0){
    uint32_t ai;
    uint32_t di;
    uint32_t si;
    std::vector<uint32_t> mi;
    std::vector<uint32_t> mi_vdc = {1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1, 1,1};
    while( offset-- && readOneInitVector(in, di, si, ai, mi) );
    int i = 0;
    // first, push a pseudo-entry for dim 0 == van der corput
    d.push_back(0);
    s.push_back(0);
    a.push_back(0);
    m.push_back(std::move(mi_vdc));
    // now, read enties from the input file
    while( i++ != nbMax && readOneInitVector(in, di, si, ai, mi) ){
        d.push_back(di);
        s.push_back(si);
        a.push_back(ai);
        if(si == 0)
        	m.push_back(std::move(mi_vdc));
        else
        	m.push_back(std::move(mi));
   }
}

inline void init_sobols(std::vector<SobolGenerator1D<uint32_t> > &sobols,
                        const std::vector<uint32_t>& d,
                        const std::vector<uint32_t>& s,
                        const std::vector<uint32_t>& a,
                        const std::vector<std::vector<uint32_t>>& m,
                        const uint32_t seed = 13374269,
				        const bool dbg_flag = false
    ) {
    sobols.resize(m.size());
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);	// uniform distribution of uint32_ts between 0 and 2^32-1
    for(int i = 0; i < m.size(); ++i) {
        sobols[i].init1D(d[i], s[i], a[i], m[i], seed, dbg_flag);
        sobols[i].owen_seeds_per_dim = unif32bits(gen);
    }
}

// extra tools

inline void loadSobolsFromFile(
        const std::string& inputFile,
        std::vector<SobolGenerator1D<uint32_t> >& sobols,
        const uint32_t seed = 13374269,
        const bool dbg_flag = false
    ) {
    std::vector<uint32_t> d;
    std::vector<uint32_t> s;
    std::vector<uint32_t> a;
    std::vector<std::vector<uint32_t>> m;
    std::ifstream tableFile(inputFile);
    if (!tableFile.is_open()) {
        std::cerr << inputFile << " cannot be read." << std::endl;
        exit(1);
    };
    load_init_table(tableFile, d, s, a, m);
    if (dbg_flag) std::cout << "loading file " << inputFile << " " << d.size() << " entries." << dbg_flag <<std::endl;
    init_sobols(sobols, d, s, a, m, seed, dbg_flag);
}	// loadSobolsFromFile

inline float getOwenPlus1D(
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
        const uint32_t dim,
        const uint32_t n,
        const uint32_t owen_tree_depth = 32,
        const bool dbg_flag = false,
        const bool owen_permut_flag = true
    ) {
    uint32_t IDcode = n;
    if (dbg_flag) std::cout << std::bitset<32>(IDcode) << " ";
    IDcode = sobols[dim].getSobolInt(n);
    if (dbg_flag) std::cout << std::bitset<32>(IDcode) << " ";
    if(owen_permut_flag) {
        IDcode = OwenScrambling(IDcode, sobols[dim].owen_seeds_per_dim, owen_tree_depth);
        if (dbg_flag) std::cout << std::bitset<32>(IDcode) << " seed:" << sobols[dim].owen_seeds_per_dim << "/" << owen_tree_depth << " ";
    }
    float res_float = ((double) IDcode / (double) UINT32SOBOLNORM);
    if (dbg_flag) std::cout << res_float << std::endl;
    return res_float;
}	// getOwenPlus1D


inline float getOwenPlus1D_with_seed(
        const std::vector<SobolGenerator1D<uint32_t> >& sobols,
        const uint32_t dim,
        const uint32_t n,
        const uint32_t owen_tree_seed,
        const uint32_t owen_tree_depth = 32,
        const bool dbg_flag = false,
        const bool owen_permut_flag = true
    ) {
    uint32_t IDcode = n;
    if (dbg_flag) std::cout << std::bitset<32>(IDcode) << " ";
    IDcode = sobols[dim].getSobolInt(n);
    if (dbg_flag) std::cout << std::bitset<32>(IDcode) << " ";
    if(owen_permut_flag) {
        IDcode = OwenScrambling(IDcode, owen_tree_seed, owen_tree_depth);
        if (dbg_flag) std::cout << std::bitset<32>(IDcode) << " seed:" << owen_tree_seed << "/" << owen_tree_depth << " ";
    }
    float res_float = ((double) IDcode / (double) UINT32SOBOLNORM);
    if (dbg_flag) std::cout << res_float << std::endl;
    return res_float;
}	// getOwenPlus1D

#endif //MYSOBOL_SOBOLGENERATOR1D_H

//std::uniform_int_distribution<uint32_t> unifFull(0);
//std::uniform_int_distribution<uint32_t> unif32bits(0, 4294967295U);	// uniform distribution of uint32_ts between 0 and 2^32-1
//std::uniform_real_distribution<double> unifreal01(0, 1.0);	// uniform distribution of reals between 0 and 1.

