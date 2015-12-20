#ifndef QSIMMATH_H
#define QSIMMATH_H

#include <gsl/gsl.h>
#include "QSimModel.h"
#include "QSimConstants.h"

namespace qsim {

    class QSimMath {
    public:

        // Make this class uninstantiable
        QSimMath() = delete;
        QSimMath(const QSimMath &math) = delete;
        QSimMath &operator=(const QSimMath &math) = delete;

        // QM Operators
        static void H    (double *f, QSimModel *model);
        static void HNorm(double *f, QSimModel *model);
        static void U    (QSimModel *model);
        static void D    (double *f, double x_range);
        static void D2   (double *f, double x_range);

        // Standard math operations
        static void add(double *f_dst, double *f_src);
        static void add(double *f, gsl_complex z);
        static void add_real(double *f, double x);
        static void add_imag(double *f, double y);
        static void sub(double *f_dst, double *f_src);
        static void multiply(double *f, gsl_complex z);
        static void multiply(double *f_dst, double *f_src);
        static void multiply_real(double *f, double x);
        static void multiply_imag(double *f, double y);
        static void multiply_add(double *f, gsl_complex z0, gsl_complex z1);
        static void multiply_add_real(double *f, double x0, double x1);
        static double max_real(const double *f);
        static double min_real(const double *f);
    };
}

// Extra macros
#ifndef MAX
#define MAX(a,b) (a>b)?a:b
#endif
#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif

#endif