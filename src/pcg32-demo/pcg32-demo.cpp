/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

/*
 * This is the original demo application from the PCG library ported to the new API
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>

#include "pcg32.h"

//--------------------------------- RNG by JC

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

int main(int argc, char** argv) {
    // Read command-line options
    int rounds = 5;

    if (argc > 1)
        rounds = atoi(argv[1]);

    pcg32 rng;

    // You should *always* seed the RNG.  The usual time to do it is the
    // point in time when you create RNG (typically at the beginning of the
    // program).
    //
    // pcg32::seed takes two 64-bit constants (the initial state, and the
    // rng sequence selector; rngs with different sequence selectors will
    // *never* have random sequences that coincide, at all)
    rng.seed(42u, 54u);

    printf("pcg32_random_r:\n"
           "      -  result:      32-bit unsigned int (uint32_t)\n"
           "      -  period:      2^64   (* 2^63 streams)\n"
           "      -  state type:  pcg32_random_t (%zu bytes)\n"
           "      -  output func: XSH-RR\n"
           "\n",
           sizeof(pcg32));

    for (int round = 1; round <= rounds; ++round) {
        printf("Round %d:\n", round);
        /* Make some 32-bit numbers */
        printf("  32bit:");
        for (int i = 0; i < 6; ++i)
            printf(" 0x%08x", rng.nextUInt());
        printf("\n");

        /* Toss some coins */
        printf("  Coins: ");
        for (int i = 0; i < 65; ++i)
            printf("%c", rng.nextUInt(2) ? 'H' : 'T');
        printf("\n");

        /* Roll some dice */
        printf("  Rolls:");
        for (int i = 0; i < 33; ++i) {
            printf(" %d", (int)rng.nextUInt(6) + 1);
        }
        printf("\n");

        /* Deal some cards */
        enum { SUITS = 4, NUMBERS = 13, CARDS = 52 };
        char cards[CARDS];

        for (int i = 0; i < CARDS; ++i)
            cards[i] = i;

        rng.shuffle(cards, cards + CARDS);

        printf("  Cards:");
        static const char number[] = {'A', '2', '3', '4', '5', '6', '7',
                                      '8', '9', 'T', 'J', 'Q', 'K'};
        static const char suit[] = {'h', 'c', 'd', 's'};
        for (int i = 0; i < CARDS; ++i) {
            printf(" %c%c", number[cards[i] / SUITS], suit[cards[i] % SUITS]);
            if ((i + 1) % 22 == 0)
                printf("\n\t");
        }
        printf("\n");

        printf("\n");
    }

    time_t start, stop;

    time(&start);
    printf("  2000000000 calls rng.seed(): \n");
    for (int i = 0; i < 2000000000; ++i) {
    	rng.seed(i, 2000000000-i);
    }
    time(&stop);
    double seconds = difftime(stop,start);
    printf("done. %f secs\n", seconds);

    RNG rngJC;

    time(&start);
    printf("  2000000000 calls rngJC.seed(): \n");
    for (int i = 0; i < 2000000000; ++i) {
    	rngJC.seed(i);
    }
    time(&stop);
    seconds = difftime(stop,start);
    printf("done. %f secs\n", seconds);

    return 0;
}
