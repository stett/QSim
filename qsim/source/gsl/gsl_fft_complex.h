#ifndef _GSL_FFT_COMPLEX_H_
#define _GSL_FFT_COMPLEX_H_

#include "gsl/gsl_complex.h"
#include <stdio.h>
#include <stdlib.h>

/**********************************************************************
 * FFT Data Types & Temporary macros
 **********************************************************************/
typedef enum
  {
    forward = -1, backward = +1,
    gsl_fft_forward = -1, gsl_fft_backward = +1      
  }
gsl_fft_direction;

/**********************************************************************
 * Function Declarations
 **********************************************************************/
int gsl_fft_complex_radix2_forward(gsl_complex_packed_array data, const size_t stride, const size_t n);
int gsl_fft_complex_radix2_backward(gsl_complex_packed_array data, const size_t stride, const size_t n);
int gsl_fft_complex_radix2_transform(gsl_complex_packed_array data, const size_t stride, const size_t n, const gsl_fft_direction sign);
int gsl_fft_complex_radix2_inverse(gsl_complex_packed_array data, const size_t stride, const size_t n);
static int fft_binary_logn (const size_t n);
static int fft_complex_bitreverse_order(gsl_complex_packed_array data, const size_t stride, const size_t n, size_t logn);

#endif