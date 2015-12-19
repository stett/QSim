#ifndef WAVEFUNCTIONPRESETS_H
#define WAVEFUNCTIONPRESETS_H

#include <gsl/gsl_complex.h>

namespace qsim {
    gsl_complex gaussian_0(double x, double x_min, double x_max);
    gsl_complex square_well_0(double x, double x_min, double x_max);
}

#endif