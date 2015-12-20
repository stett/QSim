#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>
#include "qsim/WaveFunctionPresets.h"

gsl_complex qsim::gaussian_0(double x, double x_min, double x_max) {
    double x0 = (x_min + x_max) * 0.5;
    double k0 = 20.0;
    double aa = 1.5;
    double      a = 1 / sqrt(sqrt(2 * M_PI) * aa);
    gsl_complex b = gsl_complex_exp(gsl_complex_rect(0.0, 0.0));
    //gsl_complex b = gsl_complex_exp(gsl_complex_rect(0.0, k0 * (x - x0)));
    gsl_complex c = gsl_complex_exp(gsl_complex_rect(-(x - x0) * (x - x0) / (4.0 * aa * aa), 0.0));
    gsl_complex d = gsl_complex_mul_real(b, a);
                d = gsl_complex_mul(d, c);
    return d;
}

gsl_complex qsim::square_well_0(double x, double x_min, double x_max) {
    double V_max        = 100.0;
    //double edge_width   = .45 * (x_max - x_min);
    double edge_width   = .15 * (x_max - x_min);
    return gsl_complex_rect((x < x_min + edge_width || x > x_max - edge_width) ? V_max : 0.0, 0.0);
}