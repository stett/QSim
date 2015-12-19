#include <gsl/gsl_complex.h>
#include "qsim/QSimConstants.h"
#include "qsim/QSimModel.h"
#include "qsim/QSimConstants.h"

// Coordinate conversion macros
#define QSIM_COORD_INDEX_TO_SPACE_X(n, space_min_x, space_range_x) ((double)n * space_range_x * qsim::N_INV + space_min_x)


qsim::QSimModel::QSimModel(WaveFunction psi_0, WaveFunction V_0) : psi_0(psi_0), V_0(V_0) {

    //
    // TEMP
    //
    _x_min = -1.0;
    _x_max = 10.0;
    _y_min = -1.0;
    _y_max = 2.0;
    //
    // END TEMP
    //


    // Temp vars
    gsl_complex psi_val;
    gsl_complex V_val;
    double x;

    // Loop through the points of the psi and V arrays & save the calculated values
    for (int n = 0; n < N; n ++) {

        // Find this x-coordinate
        x = QSIM_COORD_INDEX_TO_SPACE_X(n, x_min(), x_range());

        // Get the values of the psi and V functions at this point
        psi_val = psi_0(x, x_min(), x_max());
        V_val   = V_0(x, x_min(), x_max());

        // Record the values
        GSL_COMPLEX_PACKED_REAL(psi, 1, n) = GSL_REAL(psi_val);
        GSL_COMPLEX_PACKED_IMAG(psi, 1, n) = GSL_IMAG(psi_val);
        GSL_COMPLEX_PACKED_REAL(V, 1, n)   = GSL_REAL(V_val);
        GSL_COMPLEX_PACKED_IMAG(V, 1, n)   = GSL_IMAG(V_val);
    }

    // Get the initial psi_abs2
    compute_psi_abs2();
}

qsim::QSimModel::~QSimModel() {

}

void qsim::QSimModel::evolve() {

}

void qsim::QSimModel::compute_psi_abs2() {

    /*
    // Initialize the integral (sum of points) to zero
    psi_norm = 0.0;

    // Loop through each point
    double y;
    for (int n = 0; n < NN; n ++) {

        // Calculate the value at this point
        y = gsl_complex_abs(gsl_complex_rect(GSL_COMPLEX_PACKED_REAL(psi, 1, n),
                                             GSL_COMPLEX_PACKED_IMAG(psi, 1, n)));
        y = y * y;

        // Add to the sum
        psi_norm += y;

        // Save this point in the psi_abs2 array
        psi_abs2[n] = y;
    }

    // Scale the norm calculation (dx = x_range / NN)
    psi_norm *= x_range() * N_INV;
    */
}

double qsim::QSimModel::x_min() {
    return _x_min;
}

double qsim::QSimModel::x_max() {
    return _x_max;
}

double qsim::QSimModel::x_range() {
    return _x_max - _x_min;
}

double qsim::QSimModel::y_min() {
    return _y_min;
}

double qsim::QSimModel::y_max() {
    return _y_max;
}

double qsim::QSimModel::y_range() {
    return _y_max - _y_min;
}