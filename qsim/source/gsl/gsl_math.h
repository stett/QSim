#ifndef _GSL_MATH_H_
#define _GSL_MATH_H_

#include "gsl/gsl_constants.h"
#include <math.h>

/**********************************************************************
 * Function Declarations
 **********************************************************************/
double gsl_log1p (const double x);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

/**********************************************************************
 * Macros
 **********************************************************************/
#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

#endif