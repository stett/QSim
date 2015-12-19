#ifndef QSIM_H
#define QSIM_H
#include <gsl/gsl_complex.h>
#include "qsim/QSimConstants.h"

namespace qsim {
    class QSimModel {

        // Public types
    public:
        typedef gsl_complex (*WaveFunction)(double x, double x_min, double x_max);

        // Private members
    private:
        WaveFunction psi_0; // The analytical initial wavefunction, a pointer to a function of the form z=f(x)
        WaveFunction V_0;   // The analytical potential function, a pointer to a function of the form V=f(x)

        double psi[2*N];        // The wave function of a the particle
        double psi_abs2[N];     // Since abs(psi)^2 is used in multiple places, we will just always calculate it once
        double psi_norm;        // The integral of abs(psi)^2 ... should always == 1 if psi is normalized
        double V[2*N];          // The potential field
        double mass;            // The mass of the particle

        double _x_min;
        double _x_max;
        double _y_min;
        double _y_max;

        // 'Tors
    public:
        QSimModel(WaveFunction psi_0, WaveFunction V_0);
        ~QSimModel();

        // Non-copyable
    public:
        QSimModel(const QSimModel &model) = delete;
        QSimModel &operator=(const QSimModel &model) = delete;

        // Public methods
    public:
        void evolve();

        // Internal methods
    private:
        void compute_psi_abs2();

        // Getters
    public:
        double x_min();
        double x_max();
        double x_range();
    };
}

#endif