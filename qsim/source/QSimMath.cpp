#include <gsl/gsl.h>
#include "qsim/QSimMath.h"
#include "qsim/QSimModel.h"
#include "qsim/QSimConstants.h"

double qsim::QSimMath::alpha(QSimModel *model) {
    return (model->dt_size() * model->E_range()) / (2.0 * Hb);
}

gsl_complex qsim::QSimMath::lambda(QSimModel *model) {
    return gsl_complex_exp(gsl_complex_rect(0, -(model->dt_size() / HbC) * (0.5 * model->E_range() + model->E_min())));
}

int qsim::QSimMath::M(QSimModel *model) {
    double _alpha = alpha(model);
    return (int)ceil(MAX(20.0, _alpha + 11.38 * pow(_alpha, 0.32))) + 1;
}

void qsim::QSimMath::H(double *f, QSimModel *model) {

    // Make a copy of the psi array
    static double psi_0[2*N];
    for (int n = 0; n < 2*N; ++n)
        psi_0[n] = f[n];

    // Take the second derivative of the original array
    D2(f, model->x_range());

    // Loop through the points to multiply by Hb^2 / 2m and add V
    double neg_Hb2_over_2m = -HbC2 / (2.0 * model->mass());
    for (int n = 0; n < N; ++n) {

        // Calculate and save the point
        GSL_COMPLEX_PACKED_SET(f, 1, n,
            gsl_complex_add(
                gsl_complex_mul_real(GSL_COMPLEX_PACKED_GET(f, 1, n), neg_Hb2_over_2m),
                gsl_complex_mul(
                    GSL_COMPLEX_PACKED_GET(psi_0, 1, n),
                    GSL_COMPLEX_PACKED_GET(model->get_V(), 1, n)
                )
            )
        );
    }
}

void qsim::QSimMath::HNorm(double *f, QSimModel *model) {

    // Make a copy of the psi array
    static double psi_0[2*N];
    for (int n = 0; n < 2*N; ++n)
        psi_0[n] = f[n];

    // (2 / DeltaE) * (psi_0 * (-(Emin + DeltaE) / 2) + D2[psi_0] * (-Hb^2 / 2m) + psi_0 * V)
    D2(f, model->x_range());
    gsl_complex psi_times_V;
    double neg_Hb2_over_2m = -HbC2 / (2.0 * model->mass());
    double neg_Emin_minus_DeltaE_over_2 = -(model->E_min() + model->E_range() * .5);
    double two_over_DeltaE = 2./model->E_range();
    for (int n = 0; n < N; ++n) {
        psi_times_V = gsl_complex_mul(GSL_COMPLEX_PACKED_GET(psi_0, 1, n), GSL_COMPLEX_PACKED_GET(model->get_V(), 1, n));
        GSL_COMPLEX_PACKED_REAL(f, 1, n) = two_over_DeltaE * (GSL_COMPLEX_PACKED_REAL(psi_0, 1, n) * neg_Emin_minus_DeltaE_over_2 + GSL_COMPLEX_PACKED_REAL(f, 1, n) * neg_Hb2_over_2m + GSL_REAL(psi_times_V));
        GSL_COMPLEX_PACKED_IMAG(f, 1, n) = two_over_DeltaE * (GSL_COMPLEX_PACKED_IMAG(psi_0, 1, n) * neg_Emin_minus_DeltaE_over_2 + GSL_COMPLEX_PACKED_IMAG(f, 1, n) * neg_Hb2_over_2m + GSL_IMAG(psi_times_V));
    }
}

// Time evolution operator
void qsim::QSimMath::U(QSimModel *model) {

    // Make a variable to hold the coefficients for the terms in the summation
    double a_m;

    // Save the first two phi_m(x)'s
    // Zero the model's psi function
    static double phi_0[2*N];
    static double phi_1[2*N];
    static double phi_m[2*N];
    for (int n = 0; n < 2*N; ++n) {
        phi_0[n] =
        phi_1[n] =
        phi_m[n] = model->psi[n];
        model->psi[n] = 0;
    }
    HNorm(phi_1, model);
    multiply_imag(phi_1, -1.0);

    //
    double _alpha = alpha(model);

    // Begin summation loop
    for (int m = 0; m < M(model); m ++) {

        // Get phi_m
        if (m > 1) { // This is the nth time through, m > 2

            // Find the new phi_m from the previous two.
            // NOTE: phi_1 -> phi_m_minus_1
            //       phi_0 -> phi_m_minus_2
            for (int n = 0; n < 2*N; ++n)
                phi_m[n] = phi_1[n];
            HNorm(phi_m, model);
            multiply_imag(phi_m, -2.0);
            add(phi_m, phi_0);

            // Update the containers for the previous two phi_m's
            for (int n = 0; n < 2*N; ++n) {
                phi_0[n] = phi_1[n];
                phi_1[n] = phi_m[n];
            }

        } else if (m > 0) { // This is the 2nd time through
            for (int n = 0; n < 2*N; ++n)
                phi_m[n] = phi_1[n];
        }

        // Find the current expansion coefficient, a_m(t)
        a_m = (m != 0) ? 2.0 * jn(m, _alpha) : j0(_alpha);

        // Multiply this component by its coefficient & then add it to the accumulating wave packet
        multiply_real(phi_m, a_m);
        add(model->psi, phi_m);
    }

    // Multiply the entire sum by the Hnorm scaling factor, lambda(t)
    multiply(model->psi, lambda(model));
}

// First derivative operator
void qsim::QSimMath::D(double *f, double x_range) {

    // Temp vars
    int s;
    double k_over_s;

    // Calculate the component of i * k that is not dependent on the element
    k_over_s = 2. * M_PI / x_range;

    // Forward Fourier transform the function to momentum space
    gsl_fft_complex_radix2_forward(f, 1, N);

    // Multiply the first half of elements by i * k
    for (s = 0; s <= N_OVER_2; s ++) {

        // Multiply this point by i*k
        GSL_COMPLEX_PACKED_SET(f, 1, s, gsl_complex_mul_imag(GSL_COMPLEX_PACKED_GET(f, 1, s), k_over_s * s));
    }

    // Multiply the second half of elements by i * k * (s - NN) / s
    for (; s < N; s ++) {

        // Multiply this point by i*k*s/(s-NN)
        GSL_COMPLEX_PACKED_SET(f, 1, s, gsl_complex_mul_imag(GSL_COMPLEX_PACKED_GET(f, 1, s), k_over_s * (s - N)));
    }

    // Transform back to position space
    gsl_fft_complex_radix2_inverse(f, 1, N);
}

// Second derivative operator
void qsim::QSimMath::D2(double *f, double x_range) {

    // Temp vars
    int s, s_minus_N;
    double k_over_s;

    // Calculate the component of i * k that is not dependent on the element
    k_over_s = -4. * M_PI2 / (x_range * x_range);

    // Forward Fourier transform the function to momentum space
    gsl_fft_complex_radix2_forward(f, 1, N);

    // Multiply the first half of elements by -k^2
    for (s = 0; s <= N_OVER_2; s ++) {

        // Multiply this point by -k^2
        GSL_COMPLEX_PACKED_SET(f, 1, s, gsl_complex_mul_real(GSL_COMPLEX_PACKED_GET(f, 1, s), k_over_s * s * s));
    }

    // Multiply the second half of elements by -(k * (s - NN) / s)^2
    for (; s < N; s ++) {

        // Multiply this point by -(k*(s-NN)/s)^2
        s_minus_N = s - N;
        GSL_COMPLEX_PACKED_SET(f, 1, s, gsl_complex_mul_real(GSL_COMPLEX_PACKED_GET(f, 1, s), k_over_s * s_minus_N * s_minus_N));
    }

    // Transform back to position space
    gsl_fft_complex_radix2_inverse(f, 1, N);
}

/// Standard math operations
void qsim::QSimMath::add(double *f_dst, double *f_src) {
    for (int n = 0; n < N; n ++) {
        GSL_COMPLEX_PACKED_REAL(f_dst, 1, n) = GSL_COMPLEX_PACKED_REAL(f_dst, 1, n) + GSL_COMPLEX_PACKED_REAL(f_src, 1, n);
        GSL_COMPLEX_PACKED_IMAG(f_dst, 1, n) = GSL_COMPLEX_PACKED_IMAG(f_dst, 1, n) + GSL_COMPLEX_PACKED_IMAG(f_src, 1, n);
    }
}
void qsim::QSimMath::add(double *f, gsl_complex z) {
}
void qsim::QSimMath::add_real(double *f, double x) {
}
void qsim::QSimMath::add_imag(double *f, double y) {
}
void qsim::QSimMath::sub(double *f_dst, double *f_src) {
    for (int n = 0; n < N; n ++) {
        GSL_COMPLEX_PACKED_REAL(f_dst, 1, n) = GSL_COMPLEX_PACKED_REAL(f_dst, 1, n) - GSL_COMPLEX_PACKED_REAL(f_src, 1, n);
        GSL_COMPLEX_PACKED_IMAG(f_dst, 1, n) = GSL_COMPLEX_PACKED_IMAG(f_dst, 1, n) - GSL_COMPLEX_PACKED_IMAG(f_src, 1, n);
    }
}
void qsim::QSimMath::multiply(double *f, gsl_complex z) {
    for (int n = 0; n < N; n ++) {
        GSL_COMPLEX_PACKED_SET(f, 1, n, gsl_complex_mul(GSL_COMPLEX_PACKED_GET(f, 1, n), z));
    }
}
void qsim::QSimMath::multiply(double *f_dst, double *f_src) {
    for (int n = 0; n < N; n ++) {
        GSL_COMPLEX_PACKED_SET(f_dst, 1, n, gsl_complex_mul(GSL_COMPLEX_PACKED_GET(f_dst, 1, n), GSL_COMPLEX_PACKED_GET(f_src, 1, n)));
    }
}
void qsim::QSimMath::multiply_real(double *f, double x) {
    for (int n = 0; n < N; n ++) {
        GSL_COMPLEX_PACKED_REAL(f, 1, n) *= x;
        GSL_COMPLEX_PACKED_IMAG(f, 1, n) *= x;
    }
}
void qsim::QSimMath::multiply_imag(double *f, double y) {
    double x0;
    for (int n = 0; n < N; n ++) {
        x0 = GSL_COMPLEX_PACKED_REAL(f, 1, n);
        GSL_COMPLEX_PACKED_REAL(f, 1, n) = -y * GSL_COMPLEX_PACKED_IMAG(f, 1, n);
        GSL_COMPLEX_PACKED_IMAG(f, 1, n) = y * x0;
    }
}
void qsim::QSimMath::multiply_add(double *f, gsl_complex z0, gsl_complex z1) {
}
void qsim::QSimMath::multiply_add_real(double *f, double x0, double x1) {
}

// Get the minima and maxima of the real part of a function
double qsim::QSimMath::max_real(const double *f) {

    // Loop through each point and find the max
    double y;
    double y_max = 0;
    for (int n = 0; n < N; n ++) {

        // Normalize and save this point
        y = GSL_COMPLEX_PACKED_REAL(f, 1, n);
        y_max = (y > y_max) ? y : y_max;
    }
    return y_max;
}
double qsim::QSimMath::min_real(const double *f) {

    // Loop through each point and find the max
    double y;
    double y_min = 0;
    for (int n = 0; n < N; n ++) {

        // Normalize and save this point
        y = GSL_COMPLEX_PACKED_REAL(f, 1, n);
        y_min = (y < y_min) ? y : y_min;
    }
    return y_min;
}
