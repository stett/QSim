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

        // Setters
    public:
        void set_x_min(double x_min);
        void set_x_max(double x_max);
        void set_y_min(double y_min);
        void set_y_max(double y_max);

        // Getters
    public:
        double const *get_psi() const;
        double const *get_psi_abs2() const;
        double const *get_V() const;
        double x_min() const;
        double x_max() const;
        double x_range() const;
        double y_min() const;
        double y_max() const;
        double y_range() const;
    };
}

#endif