#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>
#include "qsim/WaveFunctionPresets.h"

gsl_complex qsim::Gaussian::operator()(double x) const {
    double      a = 1 / sqrt(sqrt(2 * M_PI) * alpha);
    gsl_complex b = gsl_complex_exp(gsl_complex_rect(0.0, k0 * (x - x0)));
    gsl_complex c = gsl_complex_exp(gsl_complex_rect(-(x - x0) * (x - x0) / (4.0 * alpha * alpha), 0.0));
    gsl_complex d = gsl_complex_mul_real(b, a);
    return gsl_complex_mul(d, c);
}

gsl_complex qsim::SquareWell::operator()(double x) const {
    return gsl_complex_rect((x < x_min || x > x_max) ? V_max : 0.0, 0.0);
}

gsl_complex qsim::SquareBarrier::operator()(double x) const {
    return gsl_complex_rect((x < x_min || x > x_max) ? 0.0 : V_max, 0.0);
}

/*
gsl_complex qsim::gaussian_0(double x, double x_min, double x_max) {
    double x0 = (x_min + x_max) * 0.3;
    double k0 = 10.0;
    double aa = 1.5;

    //(Sqrt[2*Pi]*a)^(-1/2)*Exp[I*k0*(x - x0)]*Exp[-(x - x0)^2/(4*a^2)]];


    double      a = 1 / sqrt(sqrt(2 * M_PI) * aa);
    gsl_complex b = gsl_complex_exp(gsl_complex_rect(0.0, k0 * (x - x0)));
    gsl_complex c = gsl_complex_exp(gsl_complex_rect(-(x - x0) * (x - x0) / (4.0 * aa * aa), 0.0));
    gsl_complex d = gsl_complex_mul_real(b, a);
                d = gsl_complex_mul(d, c);
    return d;
}

gsl_complex qsim::square_well_0(double x, double x_min, double x_max) {
    double V_max        = 16.0;
    //double edge_width   = .45 * (x_max - x_min);
    double edge_width   = .15 * (x_max - x_min);
    return gsl_complex_rect((x < x_min + edge_width || x > x_max - edge_width) ? V_max : 0.0, 0.0);
}

gsl_complex qsim::square_barrier_0(double x, double x_min, double x_max) {
    double V_max        = 4.0;
    double thickness    = 0.48 * (x_max - x_min);
    return gsl_complex_rect((x < x_min + thickness || x > x_max - thickness) ? 0.0 : V_max, 0.0);
}
*/