#ifndef WAVEFUNCTIONPRESETS_H
#define WAVEFUNCTIONPRESETS_H

#include <gsl/gsl_complex.h>

namespace qsim {

    struct Gaussian {
        double x0;// = (x_min + x_max) * 0.3;
        double k0;// = 10.0;
        double alpha;// = 1.5;
        Gaussian(double x0=0.0, double k0=10.0, double alpha=1.5) : x0(x0), k0(k0), alpha(alpha) {}
        gsl_complex operator()(double x) const;
    };

    struct SquareWell {
        double x_min;
        double x_max;
        double V_max;
        SquareWell(double x_min=0.0, double x_max=100.0, double V_max=6.0) : x_min(x_min), x_max(x_max), V_max(V_max) {}
        gsl_complex operator()(double x) const;
    };

    struct SquareBarrier {
        double x_min;
        double x_max;
        double V_max;
        SquareBarrier(double x_min=0.0, double x_max=100.0, double V_max=6.0) : x_min(x_min), x_max(x_max), V_max(V_max) {}
        gsl_complex operator()(double x) const;
    };
}

#endif