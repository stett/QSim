#ifndef QSIMCONSTANTS_H
#define QSIMCONSTANTS_H
#include <gsl/gsl_math.h>

namespace qsim {

    // The number of points of precision to keep... must be a power of 2
    const int N         = 512;
    const int N_HALF    = N/2;
    const int N_TIMES_2 = N*2;
    const int N_SQRT    = (int)sqrt(N);
    const double N_INV  = 1.0/(double)N;

    // Physical constants that are needed
    const double Hb     = 0.6582118991312933;   // [ eV fs          ]   Reduced Planck's Constant
    const double Hb2    = 0.4332429041580238;   // [(eV fs)^2       ]
    const double C      = 299.792458;           // [ nm / fs        ]   Speed of light
    const double C2     = 89875.51787368178;    // [(nm / fs)^2     ]
    const double HbC    = 197.3269631254185;    // [ eV nm          ]
    const double HbC2   = 38937.930376300275;   // [(eV nm)^2       ]
    const double EMASS  = 510998.90984764067;   // [eV fs^2 / nm^2  ]   Electron Mass

    // Constants extending those from gslp
    const double M_PI2  = M_PI * M_PI;
}

#endif