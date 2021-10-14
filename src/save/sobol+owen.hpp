// sobol+owen.hpp
//		 vo sept, 2020 essential routines for basic sobol and owen+ (no-storage permut trees)

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <bitset>
#include <list>
#include <math.h>       /* sqrt, acos etc */
#include <unistd.h>
#include <ctime>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <complex.h>
#include <iterator>
#include <algorithm>    // std::max
#include <alloca.h>
#include <limits.h>
#include <sys/time.h> // gettimeofday()
#include <random>
#include <cstdlib>

using namespace std;

//--------------------------------- constants
#define DRAND48_FLAG // drand48()
#undef DRAND48_FLAG // drand48()  -- don't use it -- it's not good

#define MX_SIZE 32
#define N_SIGNIFICANT_BITS 32

#define N_SOBOL_INIT_TAB_ENTRIES    100000 	// KJ table: 21201 but we can potentially handle more 

#define UINT32NORM 4294967296U
#define UINT64NORM 18446744073709551615ULL

//--------------------------------- global variables

uint32_t sobol_dj[N_SOBOL_INIT_TAB_ENTRIES];
uint32_t sobol_sj[N_SOBOL_INIT_TAB_ENTRIES];
uint32_t sobol_aj[N_SOBOL_INIT_TAB_ENTRIES];
uint32_t sobol_mk[N_SOBOL_INIT_TAB_ENTRIES][CHAR_BIT * sizeof(uint32_t)];
uint32_t powers_of_two[N_SIGNIFICANT_BITS];

//--------------------------------- structures
typedef unsigned char uchar;

struct t_double_point2D {
	double x, y;
};
struct t_uint_point2D {
	uint32_t x, y;
};
struct t_int_point2D {
	int32_t x, y;
};

struct t_sampler_env {
	uint32_t nDims;
	uint32_t *sobolIndices;
	uint32_t *randomSeed_Binary; //one random seed per dim
	uint32_t rng_state;
};


//--------------------------------- routines

//--------------------------------- core routines
//--------------------------------- Sobol-related routines
// ReverseBits32() -- to reverse bits [0..31]
// ReverseBits64() -- to reverse bits [0..63]

inline uint32_t ReverseBits32(uint32_t n) {
	n = (n << 16) | (n >> 16);
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
	return n;
}

inline uint64_t ReverseBits64(uint64_t n) {
	n = (n << 32) | (n >> 32);
	n = ((n & 0x0000ffff0000ffffULL) << 16) | ((n & 0xffff0000ffff0000ULL) >> 16);
	n = ((n & 0x00ff00ff00ff00ffULL) << 8) | ((n & 0xff00ff00ff00ff00ULL) >> 8);
	n = ((n & 0x0f0f0f0f0f0f0f0fULL) << 4) | ((n & 0xf0f0f0f0f0f0f0f0ULL) >> 4);
	n = ((n & 0x3333333333333333ULL) << 2) | ((n & 0xccccccccccccccccULL) >> 2);
	n = ((n & 0x5555555555555555ULL) << 1) | ((n & 0xaaaaaaaaaaaaaaaaULL) >> 1);
	return n;
}

int fileExists(char *fname) {
	if (access(fname, F_OK) != -1) {
		return 1;
	} else {
		cout << fname << " file doesn't exist\n";
		return 0;
	}
}

//--------------------------------- end of general-purpose routines

//--------------------------------- TRNG
template<unsigned int rng_a, unsigned int rng_c, unsigned int rng_m>
struct TRNG {
	TRNG() :
			mstate(1u), mseed(1u) {
	}

	TRNG(const unsigned int _seed) :
			mstate(_seed), mseed(_seed) {
	}

	void seed(const unsigned int _seed) {
		//~ mseed= _seed;
		mseed = hash(_seed);
		mstate = mseed;
	}

	void index(const int _index) {
		// advance
		// cf http://www.pcg-random.org, c++ implementation
		unsigned int cur_mul = rng_a;
		unsigned int cur_add = rng_c;
		unsigned int acc_mul = 1u;
		unsigned int acc_add = 0u;
		unsigned int delta = _index;
		while (delta > 0u) {
			if (delta & 1u) {
				acc_mul = acc_mul * cur_mul;
				acc_add = acc_add * cur_mul + cur_add;
			}

			//~ cur_add= (cur_mul+1u) * cur_add;
			cur_add = cur_mul * cur_add + cur_add;
			cur_mul = cur_mul * cur_mul;
			delta = delta >> 1u;
		}

		// advance current state
		// mstate= acc_mul * mstate + acc_add;

		// advance to sample index
		mstate = acc_mul * mseed + acc_add;
	}

	unsigned int operator()() {
		return sample();
	}

	unsigned int sample() {
		mstate = (mstate * rng_a + rng_c) % rng_m;
		return mstate;
	}

	double sample_double() {
		return sample() / double(rng_m);
	}

	unsigned int sample_range(const unsigned int range) {
		// Efficiently Generating a Number in a Range
		// cf http://www.pcg-random.org/posts/bounded-rands.html
		unsigned int divisor = rng_m / range;
		if (divisor == 0)
			return 0;

		while (true) {
			unsigned int x = sample() / divisor;
			if (x < range)
				return x;
		}
	}

protected:
	unsigned int mstate;
	unsigned int mseed;

	// cf http://www.burtleburtle.net/bob/hash/integer.html
	uint32_t hash(uint32_t a) {
		a = (a + 0x7ed55d16) + (a << 12);
		a = (a ^ 0xc761c23c) ^ (a >> 19);
		a = (a + 0x165667b1) + (a << 5);
		a = (a + 0xd3a2646c) ^ (a << 9);
		a = (a + 0xfd7046c5) + (a << 3);
		a = (a ^ 0xb55a4f09) ^ (a >> 16);
		return a;
	}
};

typedef TRNG<1103515245u, 12345u, 1u << 31> RNG;
//--------------------------------- end of RNG

/* `sobol_mk` stores the direction numbers for all the dimensions. For
 * example, the direction numbers for dimension `d` are stored in the
 * memory space:
 *
 *      sobol_mk[B*d] ... sobol_mk[B*d + B-1]
 *
 * where B is the number of bits of the index.
 *
 * For 32 bits indexing the number of direction numbers per dimension
 * is 32.
 *
 * `sobol_aj' stores the primitive polynomials.
 * `sobol_sj` stores the dimension of the primitive polynomials `sobol_aj`.
 */
//uint32_t sobol_mk[21201][CHAR_BIT*sizeof(uint32_t)];
//uint32_t sobol_aj[21201];
//uint32_t sobol_sj[21201];
/* Generate the direction numbers for dimension `dim` using the provided
 * polynomial `aj` of degree `sj` and the s initial direction numbers stored
 * in the first `sj` elements of vector `mk`.
 *
 * We assume that mk is already allocated to the number of bits composing
 * the integer type.
 */
inline
void generate_mk(uint32_t aj, uint32_t sj, uint32_t* mk) {
	/* For each direction number, we apply the recurrence formula:
	 *    m_k = aj_1 2 m_{k-1} xor ... xor aj_{sj} 2^{sj} m_{k-sj}
	 *                                 xor mk_{k-sj}
	 */
	const uint32_t size = CHAR_BIT * sizeof(uint32_t);
	for (int k = sj; k < size; ++k) {
		mk[k] = 0;
		for (int i = 1; i < sj; ++i) {
			// `akj` stores aj_k, note that the polynomial rep is reversed
			// `pw2` stores 2^k
			const uint32_t akj = (aj >> (sj - 1 - i)) & 1;
			const uint32_t pw2 = (1 << i);
			mk[k] ^= akj * pw2 * mk[k - i];
		}
		mk[k] ^= (1 << sj) * mk[k - sj];
		mk[k] ^= mk[k - sj];

	}
}


/* Load the direction number, primite polynomials and associated dimensions
 * from a file. Depending on the file, we can load either the Joe and Kuo
 * file or a binary one.
 */
inline
uint32_t load_mk(const string filename, const bool dbg_flag) {
	int size = 0, dim_from = 1;
	if (dbg_flag) cout << "Loading J&K file " << filename << " ... dim_from=" << dim_from << endl;
	std::ifstream file(filename, std::ifstream::binary);
	if (!file.is_open()) {
		std::cout << "load_mk: " << filename << " not found." << "\n";
		exit(1);
	};
	/* Check if the file is supported
	 * We support two kind of file format: the one of Joe and Kuo that
	 * starts with "d\ts\ta" or one own binary format that starts with
	 * "SOBOL".
	 */
	char c = file.get();
	if(c == 'd') {	// JKheader_flag : skip header, otherwise continue loading file in Joe & Kuo format
		file.ignore(256, '\n');
	} else {
		file.putback(c);
	}
	// Load all the mk, aj and sj
	// We do not store more than 1110 dimension
	uint32_t dim = dim_from;
	while (file.good() && dim < N_SOBOL_INIT_TAB_ENTRIES) {
		// Parse the polynomial
		uint32_t d, sj, aj;
		file >> d >> sj >> aj;
		sobol_aj[dim] = aj;
		sobol_sj[dim] = sj;
		sobol_dj[dim] = d;
		// Parse the direction numbers
		for (int i = 0; i < sj; ++i) {
			file >> sobol_mk[dim][i];
		}

		// Fill the remaining direction numbers
		generate_mk(aj, sj, sobol_mk[dim]);
		dim++;
	}
	return dim-2;
}	// load_mk

/* Save the direction number, primite polynomials and associated dimensions
 * from a file.
 */
inline
void save_mk(const char* filename) {
	cout << "save mk data into " << filename << endl;
	std::ofstream file(filename, std::ifstream::trunc);
	// Print header
	file << "d\ts\ta\tm_i" << std::endl;

	// Print content
	for (int d = 0; d < N_SOBOL_INIT_TAB_ENTRIES; ++d) {
		file << sobol_dj[d] << "\t" << sobol_sj[d] << "\t" << sobol_aj[d] << "\t";
		for (int k = 0; k < sobol_sj[d]; ++k) {
			file << sobol_mk[d][k] << " ";
		}
		file << std::endl;
	}
}

/* Generate the binary representation of the `i`-th point of the Sobol
 * sequence using the `mk` direction number.
 */
inline uint32_t sobol_binary(uint32_t i, const uint32_t* mk) {
	// `x` is the result binary representation of a float
	uint32_t x = 0;
	// For the k-th bit of the index, if this bit is one, accumulate
	// the k-th direction number modulo 2 (using a xor).
	for (int k = 0; k < CHAR_BIT * sizeof(uint32_t); ++k) {
		x <<= 1;
		if (i & 1) {
			x ^= *(mk++);
		} else {
			mk++;
		}
		i >>= 1;
	}
	return x;
}

/* Generate a floating point value from 'x' */
inline
float binary_to_float(uint32_t x) {
	// `r` is the floatting point result
	float r = 0.0f;

	// For the k-th bit in x, accumulate 2^-k in `r`.
	float b = 0.5f; // 1 / 2^k
	for (int k = CHAR_BIT * sizeof(uint32_t) - 1; k >= 0; --k) {
		// Add the base if the float representation is 1 at the
		// selector
		r += ((x >> k) & 1) ? b : 0.0f;

		// Scale the base
		b *= 0.5f;
	}
	return r;
}

/* Compute the `i`-th element of the Sobol sequence using the `mk`
 * direction numbers.
 */
inline
float sobol(uint32_t i, const uint32_t* mk) {
	// Evaluate Sobol in binary form
	const uint32_t x = sobol_binary(i, mk);

	// Convert binary form to floats
	const float v = binary_to_float(x);
	return v;
}

/* Compute the `i`-th element of the Sobol sequence for the `d`-th
 * dimension. This function uses the direction numbers stored in the
 * `sobol_mk` array.
 */
inline uint32_t get_1pt_sobol1D(const uint32_t dim, const uint32_t n) {
	if (dim == 0)                // dim 0: van der Corput
		return ReverseBits32(n);
	else
		return sobol_binary(n, sobol_mk[dim]);	// 1 -> mk[1]] etc.
}

void read_mk(std::string fname_str, uint32_t *sobol_mk) {
	std::ifstream infile;
	std::cout << "reading " << fname_str << "\n";
	infile.open(fname_str, std::ios::in);
	if (!infile.is_open()) {
		std::cout << fname_str << " not found." << "\n";
		exit(1);
	};
	for (int i = 0; i < MX_SIZE; i++) {
		uint32_t current_val;
		infile >> current_val;
		sobol_mk[i] = current_val;
	}
	infile.close();
}

uint32_t getIDcodePermutedBinary(
		t_sampler_env *sampler_env,
		const uint32_t IDcode,
		const uint32_t dim,
		const uint32_t owen_tree_depth) {
	std::bitset<N_SIGNIFICANT_BITS> digitsBase2 = IDcode;
	std::bitset<N_SIGNIFICANT_BITS> new_digitsBase2 = 0; // attention: numbering in the bitset is "reversed" : 0 corresponds to the least significant bit, 31 - to the most signifcant bit
	uint32_t permut;
	uint32_t thisDigitPermut;
	RNG rng;
	for (int idigit = 0; idigit < owen_tree_depth; idigit++) {
		uint32_t indPermut = powers_of_two[idigit] - 1 + (digitsBase2 >> (N_SIGNIFICANT_BITS - idigit)).to_ulong(); // address in the tree
		uint32_t seed = sampler_env->randomSeed_Binary[dim] + (indPermut);
#ifdef DRAND48_FLAG
		srand48(seed);
		thisDigitPermut = round(drand48());
#else
		rng.seed(seed);
		thisDigitPermut = rng.sample_range(2);
#endif
		new_digitsBase2[(N_SIGNIFICANT_BITS - 1) - idigit] = (thisDigitPermut ^ digitsBase2[(N_SIGNIFICANT_BITS - 1) - idigit]);
	}
	return new_digitsBase2.to_ulong();
}    // getIDcodePermutedBinary

//--------------------------------- end of core routines

//-------------------------------- init & free routines

// sampler_init() creates all needed internal structures
t_sampler_env *sampler_init(const string dir_vectors_fname, uint32_t global_seed, const bool dbg_flag) {
	struct timeval tp;

	uint32_t nDims = load_mk(dir_vectors_fname, dbg_flag);

	cout << dir_vectors_fname << " : " << nDims << " entries read" << endl;

	struct t_sampler_env *sampler_env = (t_sampler_env*) malloc(sizeof(t_sampler_env));
	sampler_env->nDims = nDims;
	std::cout << "sampler_init: nDims=" << nDims << std::endl;
	sampler_env->randomSeed_Binary = (uint32_t*) malloc(sizeof(uint32_t) * nDims);
	for (int idim = 1; idim <= nDims; idim++) {    // skip entry [0] because it's van der Corput
		gettimeofday(&tp, NULL);
		uint32_t timeStamp = time(NULL) + tp.tv_usec;
//		if (dbg_flag) cout << idim << " timeStamp: " << timeStamp << endl;
		uint32_t idim32 = ((uint32_t) idim << 16) + ((uint32_t) idim << 8) + (uint32_t) idim;
		if (global_seed == 0) {
			global_seed = tp.tv_usec;
			srand48((idim+1)*137 + 234567);
			uint32_t seed = round(209567843*drand48());
			sampler_env->randomSeed_Binary[idim] = (idim+1)*65497 ^ seed;
		} else {
			srand48((idim+1)*5753 + 856479);
			uint32_t seed = round(168543723*drand48());
			sampler_env->randomSeed_Binary[idim] = global_seed ^ seed;
		}
		if (dbg_flag)
			std::cout << "sampler_env->randomSeed_Binary[" << idim << "] = " <<   sampler_env->randomSeed_Binary[idim] << endl;
	}
	if (dbg_flag) {
		for (int idim = 1; idim <= nDims; idim++) {    // skip entry [0] because it's van der Corput
			cout << idim << " " << dir_vectors_fname << " : d/s/a -> \t" << sobol_dj[ idim ] << "   \t" << sobol_sj[ idim ] << "   \t" << sobol_aj[ idim ] << "   \t" ;
			for (int k = 0; k < sobol_sj[ idim ]; ++k) {
				cout  << " " << sobol_mk[ idim ][ k ];
			}
			cout << endl;
		}
		for (int idim = 1; idim <= nDims; idim++) {    // skip entry [0] because it's van der Corput
			cout << idim << " " << dir_vectors_fname << " : d/s/a -> \t" << sobol_dj[ idim ] << "   \t" << sobol_sj[ idim ] << "   \t" << sobol_aj[ idim ] << "   \t" ;
			for (int k = 0; k < 32; ++k) {
				cout  << " " << sobol_mk[ idim ][ k ];
			}
			cout << endl;
		}
	}
	sampler_env->rng_state = 0;
	powers_of_two[0] = 1;
	for (int ipow = 1; ipow < N_SIGNIFICANT_BITS; ipow++)
		powers_of_two[ipow] = 2 * powers_of_two[ipow - 1];
	return sampler_env;
} // sampler_init

// sampler_init_seeds() initiates sampler_env->randomSeed_Binary[idim]
inline void sampler_init_seeds_one_entry(
		const t_sampler_env *sampler_env,
		uint32_t global_seed,
		const int idim,
		const bool dbg_flag) {
	sampler_env->randomSeed_Binary[idim] = global_seed;
//	if (dbg_flag) cout << "sampler_env->randomSeed_Binary[" << idim << "] = " <<   sampler_env->randomSeed_Binary[idim] << endl;
} // sampler_init_seeds_one_entry

//--------------------------------- end of Sobol-related routines
