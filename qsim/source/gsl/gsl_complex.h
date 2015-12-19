#ifndef _GSL_COMPLEX_H_
#define _GSL_COMPLEX_H_

#include <math.h>
#include "gsl/gsl_math.h"

/**********************************************************************
 * Complex Data Types & Macros
 **********************************************************************/

/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

typedef const double *       gsl_const_complex_packed ;
typedef const float *        gsl_const_complex_packed_float  ;
typedef const long double *  gsl_const_complex_packed_long_double ;

/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;

typedef const double *       gsl_const_complex_packed_array ;
typedef const float *        gsl_const_complex_packed_array_float  ;
typedef const long double *  gsl_const_complex_packed_array_long_double ;

/* Yes... this seems weird. Trust us. The point is just that
   sometimes you want to make it obvious that something is
   an output value. The fact that it lacks a 'const' may not
   be enough of a clue for people in some contexts.
 */
typedef double *       gsl_complex_packed_ptr ;
typedef float *        gsl_complex_packed_float_ptr  ;
typedef long double *  gsl_complex_packed_long_double_ptr ;

typedef const double *       gsl_const_complex_packed_ptr ;
typedef const float *        gsl_const_complex_packed_float_ptr  ;
typedef const long double *  gsl_const_complex_packed_long_double_ptr ;

typedef struct {
    long double dat[2];
} gsl_complex_long_double;

typedef struct {
    double dat[2];
} gsl_complex;

typedef struct {
    float dat[2];
} gsl_complex_float;

#define GSL_REAL(z)             ((z).dat[0])
#define GSL_IMAG(z)             ((z).dat[1])
#define GSL_COMPLEX_P(zp)       ((zp)->dat)
#define GSL_COMPLEX_P_REAL(zp)  ((zp)->dat[0])
#define GSL_COMPLEX_P_IMAG(zp)  ((zp)->dat[1])
#define GSL_COMPLEX_EQ(z1,z2)   (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))
/// I made the following two macros to replace a GSL function that I couldn't seem to find in the source
///  - steven
#define GSL_COMPLEX_PACKED_REAL(zp,stride,i) *((zp)+2*(stride)*(i))
#define GSL_COMPLEX_PACKED_IMAG(zp,stride,i) *((zp)+2*(stride)*(i)+1)
#define GSL_COMPLEX_PACKED_GET(zp,stride,i) gsl_complex_rect(GSL_COMPLEX_PACKED_REAL(zp,stride,i), GSL_COMPLEX_PACKED_IMAG(zp,stride,i))
#define GSL_COMPLEX_PACKED_SET(zp,stride,i,z) GSL_COMPLEX_PACKED_REAL(zp,stride,i)=GSL_REAL(z);GSL_COMPLEX_PACKED_IMAG(zp,stride,i)=GSL_IMAG(z)
/// I redefined them since the only times that they are used is with stride=1
//#define GSL_COMPLEX_PACKED_GET(zp,i) gsl_complex_rect(GSL_COMPLEX_PACKED_REAL(zp,i), GSL_COMPLEX_PACKED_IMAG(zp,i))
//#define GSL_COMPLEX_PACKED_REAL(zp,i) *((zp)+2*(i))
//#define GSL_COMPLEX_PACKED_IMAG(zp,i) *((zp)+2*(i)+1)

#define GSL_SET_COMPLEX(zp,x,y) do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#define GSL_SET_REAL(zp,x)      do {(zp)->dat[0]=(x);} while(0)
#define GSL_SET_IMAG(zp,y)      do {(zp)->dat[1]=(y);} while(0)

#define GSL_SET_COMPLEX_PACKED(zp,n,x,y) do {*((zp)+2*(n))=(x); *((zp)+(2*(n)+1))=(y);} while(0)

/**********************************************************************
 * Function Declarations
 **********************************************************************/
gsl_complex gsl_complex_rect (double x, double y);
gsl_complex gsl_complex_polar (double r, double theta);
double gsl_complex_arg (gsl_complex z);
double gsl_complex_abs (gsl_complex z);
double gsl_complex_abs2 (gsl_complex z);
double gsl_complex_logabs (gsl_complex z);
gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_add_real (gsl_complex a, double x);
gsl_complex gsl_complex_add_imag (gsl_complex a, double y);
gsl_complex gsl_complex_sub (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_sub_real (gsl_complex a, double x);
gsl_complex gsl_complex_sub_imag (gsl_complex a, double y);
gsl_complex gsl_complex_mul (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_mul_real (gsl_complex a, double x);
gsl_complex gsl_complex_mul_imag (gsl_complex a, double y);
gsl_complex gsl_complex_div (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_div_real (gsl_complex a, double x);
gsl_complex gsl_complex_div_imag (gsl_complex a, double y);
gsl_complex gsl_complex_conjugate (gsl_complex a);
gsl_complex gsl_complex_negative (gsl_complex a);
gsl_complex gsl_complex_inverse (gsl_complex a);
gsl_complex gsl_complex_sqrt (gsl_complex a);
gsl_complex gsl_complex_sqrt_real (double x);
gsl_complex gsl_complex_exp (gsl_complex a);
gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_pow_real (gsl_complex a, double b);
gsl_complex gsl_complex_log (gsl_complex a);
gsl_complex gsl_complex_log10 (gsl_complex a);
gsl_complex gsl_complex_log_b (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_sin (gsl_complex a);
gsl_complex gsl_complex_cos (gsl_complex a);
gsl_complex gsl_complex_tan (gsl_complex a);
gsl_complex gsl_complex_sec (gsl_complex a);
gsl_complex gsl_complex_csc (gsl_complex a);
gsl_complex gsl_complex_cot (gsl_complex a);
gsl_complex gsl_complex_arcsin (gsl_complex a);
gsl_complex gsl_complex_arcsin_real (double a);
gsl_complex gsl_complex_arccos (gsl_complex a);
gsl_complex gsl_complex_arccos_real (double a);
gsl_complex gsl_complex_arctan (gsl_complex a);
gsl_complex gsl_complex_arcsec (gsl_complex a);
gsl_complex gsl_complex_arccsc (gsl_complex a);
gsl_complex gsl_complex_arccsc_real (double a);
gsl_complex gsl_complex_arccot (gsl_complex a);
gsl_complex gsl_complex_sinh (gsl_complex a);
gsl_complex gsl_complex_cosh (gsl_complex a);
gsl_complex gsl_complex_tanh (gsl_complex a);
gsl_complex gsl_complex_sech (gsl_complex a);
gsl_complex gsl_complex_csch (gsl_complex a);
gsl_complex gsl_complex_coth (gsl_complex a);
gsl_complex gsl_complex_arcsinh (gsl_complex a);
gsl_complex gsl_complex_arccosh (gsl_complex a);
gsl_complex gsl_complex_arccosh_real (double a);
gsl_complex gsl_complex_arctanh (gsl_complex a);
gsl_complex gsl_complex_arctanh_real (double a);
gsl_complex gsl_complex_arcsech (gsl_complex a);
gsl_complex gsl_complex_arccsch (gsl_complex a);
gsl_complex gsl_complex_arccoth (gsl_complex a);

#endif