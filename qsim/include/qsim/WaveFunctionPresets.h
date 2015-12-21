#ifndef WAVEFUNCTIONPRESETS_H
#define WAVEFUNCTIONPRESETS_H

#include <gsl/gsl_complex.h>

namespace qsim {
    //gsl_complex gaussian_0(double x, double x_min, double x_max);
    //gsl_complex square_well_0(double x, double x_min, double x_max);
    //gsl_complex square_barrier_0(double x, double x_min, double x_max);

    struct Gaussian {
        double x0;// = (x_min + x_max) * 0.3;
        double k0;// = 10.0;
        double alpha;// = 1.5;
        Gaussian(double x0, double k0, double alpha) : x0(x0), k0(k0), alpha(alpha) {}
        gsl_complex operator()(double x) const;
    };

    struct SquareWell {
        double x_min;
        double x_max;
        double V_max;
        SquareWell(double x_min, double x_max, double V_max) : x_min(x_min), x_max(x_max), V_max(V_max) {}
        gsl_complex operator()(double x) const;
    };

    struct SquareBarrier {
        double x_min;
        double x_max;
        double V_max;
        SquareBarrier(double x_min, double x_max, double V_max) : x_min(x_min), x_max(x_max), V_max(V_max) {}
        gsl_complex operator()(double x) const;
    };
}

#endif