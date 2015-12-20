#include <stdexcept>
#include <gsl/gsl_complex.h>
#include "qsim/QSimConstants.h"
#include "qsim/QSimModel.h"
#include "qsim/QSimMath.h"
#include "qsim/QSimConstants.h"
#include "qsim/QSimCoordinates.h"


qsim::QSimModel::QSimModel(WaveFunction psi_0, WaveFunction V_0) : psi_0(psi_0), V_0(V_0) {

    //
    // TEMP
    //
    _x_min = -1.0;
    _x_max = 10.0;
    _mass = EMASS;
    _dt_size = 0.0061;
    _dt_iterations = 1;
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

qsim::QSimModel::~QSimModel() {}

void qsim::QSimModel::evolve() {

    // Update the wave form
    for (int i = 0; i < _dt_iterations; ++i)
        QSimMath::U(this);

    // Update the squared absolute value of the waveform
    compute_psi_abs2();
}

void qsim::QSimModel::compute_psi_abs2() {

    // Initialize the integral (sum of points) to zero
    psi_norm = 0.0;

    // Loop through each point
    double y;
    for (int n = 0; n < N; n ++) {

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
}

void qsim::QSimModel::set_dt_size(double dt_size) {
    if (dt_size <= 0.0f)
        throw std::invalid_argument("dt_size");
    _dt_size = dt_size;
}
void qsim::QSimModel::set_dt_iterations(double dt_iterations) {
    if (dt_iterations <= 0.0f)
        throw std::invalid_argument("dt_iterations");
}
void qsim::QSimModel::set_x_min(double x_min) { _x_min = x_min; }
void qsim::QSimModel::set_x_max(double x_max) { _x_max = x_max; }
double const *qsim::QSimModel::get_psi() const { return psi; }
double const *qsim::QSimModel::get_psi_abs2() const { return psi_abs2; }
double const *qsim::QSimModel::get_V() const { return V; }
double qsim::QSimModel::mass() const { return _mass; }
double qsim::QSimModel::dt_size() const { return _dt_size; }
double qsim::QSimModel::dt_iterations() const { return _dt_iterations; }
double qsim::QSimModel::x_min() const { return _x_min; }
double qsim::QSimModel::x_max() const { return _x_max; }
double qsim::QSimModel::x_range() const { return _x_max - _x_min; }

double qsim::QSimModel::K_max() const {
    double _K_max;
    _K_max = M_PI * HbC * N / _x_range;
    _K_max = _K_max * _K_max / (2 * _mass);
    return _K_max;
}

double qsim::QSimModel::E_min() const {
    double V_min = QSimMath::min_real(V);
    double _E_min = (0 < V_min) ? 0 : V_min;
    return _E_min;
}

double qsim::QSimModel::E_max() const {
    double V_max = QSimMath::max_real(V);
    double _E_max = K_max() + V_max;
    return _E_max;
}

double qsim::QSimModel::E_range() const {
    double _E_range = E_max() - E_min();
    return _E_range;
}
