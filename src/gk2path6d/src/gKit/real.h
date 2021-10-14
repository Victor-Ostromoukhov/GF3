
#ifndef _REAL_H
#define _REAL_H

#include <cfloat>

#ifdef FLOAT32
#warning using float32

typedef float Float;
#define FLOAT_MIN FLT_MIN
#define FLOAT_MAX FLT_MAX

#else

typedef double Float;
#define FLOAT_MIN DBL_MIN
#define FLOAT_MAX DBL_MAX
#endif

#endif
