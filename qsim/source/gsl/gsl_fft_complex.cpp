#include "gsl/gsl_fft_complex.h"

/**********************************************************************
 * Radix-2 (N = n^2)
 **********************************************************************/
int gsl_fft_complex_radix2_forward(gsl_complex_packed_array data,
                                    const size_t stride,
                                    const size_t n)
{
  gsl_fft_direction sign = gsl_fft_forward;
  int status = gsl_fft_complex_radix2_transform(data, stride, n, sign);
  return status;
}

int gsl_fft_complex_radix2_backward(gsl_complex_packed_array data,
                                    const size_t stride,
                                    const size_t n)
{
  gsl_fft_direction sign = gsl_fft_backward;
  int status = gsl_fft_complex_radix2_transform(data, stride, n, sign);
  return status;
}

int gsl_fft_complex_radix2_transform(gsl_complex_packed_array data,
                                     const size_t stride, 
                                     const size_t n,
                                     const gsl_fft_direction sign)
{
  int result ;
  size_t dual;
  size_t bit; 
  size_t logn = 0;
  int status;

  if (n == 1) /* identity operation */
    {
      return 0 ;
    }

  /* make sure that n is a power of 2 */

  result = fft_binary_logn(n) ;

  if (result == -1) 
    {
      // TODO: Print message:
      // "n is not a power of 2";
      abort();
    } 
  else 
    {
      logn = result ;
    }

  /* bit reverse the ordering of input data for decimation in time algorithm */
  
  status = fft_complex_bitreverse_order(data, stride, n, logn);

  /* apply fft recursion */

  dual = 1;

  for (bit = 0; bit < logn; bit++)
    {
      double w_real = 1.0;
      double w_imag = 0.0;

      const double theta = 2.0 * ((int) sign) * M_PI / (2.0 * (double) dual);

      const double s = sin (theta);
      const double t = sin (theta / 2.0);
      const double s2 = 2.0 * t * t;

      size_t a, b;

      /* a = 0 */

      for (b = 0; b < n; b += 2 * dual)
        {
          const size_t i = b ;
          const size_t j = b + dual;

          const double z1_real = GSL_COMPLEX_PACKED_REAL(data,stride,j);
          const double z1_imag = GSL_COMPLEX_PACKED_IMAG(data,stride,j);

          const double wd_real = z1_real ;
          const double wd_imag = z1_imag ;
          
          GSL_COMPLEX_PACKED_REAL(data,stride,j) = GSL_COMPLEX_PACKED_REAL(data,stride,i) - wd_real;
          GSL_COMPLEX_PACKED_IMAG(data,stride,j) = GSL_COMPLEX_PACKED_IMAG(data,stride,i) - wd_imag;
          GSL_COMPLEX_PACKED_REAL(data,stride,i) += wd_real;
          GSL_COMPLEX_PACKED_IMAG(data,stride,i) += wd_imag;
        }
      
      /* a = 1 .. (dual-1) */

      for (a = 1; a < dual; a++)
        {

          /* trignometric recurrence for w-> exp(i theta) w */

          {
            const double tmp_real = w_real - s * w_imag - s2 * w_real;
            const double tmp_imag = w_imag + s * w_real - s2 * w_imag;
            w_real = tmp_real;
            w_imag = tmp_imag;
          }

          for (b = 0; b < n; b += 2 * dual)
            {
              const size_t i = b + a;
              const size_t j = b + a + dual;

              const double z1_real = GSL_COMPLEX_PACKED_REAL(data,stride,j);
              const double z1_imag = GSL_COMPLEX_PACKED_IMAG(data,stride,j);
              
              const double wd_real = w_real * z1_real - w_imag * z1_imag;
              const double wd_imag = w_real * z1_imag + w_imag * z1_real;

              GSL_COMPLEX_PACKED_REAL(data,stride,j) = GSL_COMPLEX_PACKED_REAL(data,stride,i) - wd_real;
              GSL_COMPLEX_PACKED_IMAG(data,stride,j) = GSL_COMPLEX_PACKED_IMAG(data,stride,i) - wd_imag;
              GSL_COMPLEX_PACKED_REAL(data,stride,i) += wd_real;
              GSL_COMPLEX_PACKED_IMAG(data,stride,i) += wd_imag;
            }
        }
      dual *= 2;
    }

  return 0;
}

int gsl_fft_complex_radix2_inverse(gsl_complex_packed_array data,
                                   const size_t stride,
                                   const size_t n)
{
  gsl_fft_direction sign = gsl_fft_backward;
  int status = gsl_fft_complex_radix2_transform(data, stride, n, sign);

  if (status)
    {
      return status;
    }

  /* normalize inverse fft with 1/n */

  {
    const double norm = 1.0 / n;
    size_t i;
    for (i = 0; i < n; i++)
      {
        GSL_COMPLEX_PACKED_REAL(data,stride,i) *= norm;
        GSL_COMPLEX_PACKED_IMAG(data,stride,i) *= norm;
      }
  }

  return status;
}

static int fft_binary_logn (const size_t n)
{
  size_t ntest ;
  size_t binary_logn = 0 ;
  size_t k = 1;

  while (k < n)
    {
      k *= 2;
      binary_logn++;
    }

  ntest = (1 << binary_logn) ;

  if (n != ntest )       
    {
      return -1 ; /* n is not a power of 2 */
    } 

  return binary_logn;
}

static int fft_complex_bitreverse_order(gsl_complex_packed_array data,
                                        const size_t stride,
                                        const size_t n,
                                        size_t logn)
{
  /* This is the Goldrader bit-reversal algorithm */

  size_t i;
  size_t j = 0;

  logn = 0 ; /* not needed for this algorithm */

  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;

      if (i < j)
        {
          const double tmp_real = GSL_COMPLEX_PACKED_REAL(data,stride,i);
          const double tmp_imag = GSL_COMPLEX_PACKED_IMAG(data,stride,i);

          GSL_COMPLEX_PACKED_REAL(data,stride,i) = GSL_COMPLEX_PACKED_REAL(data,stride,j);
          GSL_COMPLEX_PACKED_IMAG(data,stride,i) = GSL_COMPLEX_PACKED_IMAG(data,stride,j);
          GSL_COMPLEX_PACKED_REAL(data,stride,j) = tmp_real;
          GSL_COMPLEX_PACKED_IMAG(data,stride,j) = tmp_imag;
        }

      while (k <= j) 
        {
          j = j - k ;
          k = k / 2 ;
        }

      j += k ;
    }

  return 0;
}