#ifndef QSIMMODEL_H
#define QSIMMODEL_H
#include <gsl/gsl_complex.h>
#include "qsim/QSimConstants.h"

namespace qsim {
    class QSimModel {

        // Make friends with the QSim math class
        friend class QSimMath;

        // Public types
    public:
        typedef gsl_complex (*WaveFunction)(double x, double x_min, double x_max);

        // Private members
    private:
        WaveFunction psi_0;     // The analytical initial wavefunction, a pointer to a function of the form z=f(x)
        WaveFunction V_0;       // The analytical potential function, a pointer to a function of the form V=f(x)
        double psi[2*N];        // The wave function of a the particle
        double psi_abs2[N];     // Since abs(psi)^2 is used in multiple places, we will just always calculate it once
        double psi_norm;        // The integral of abs(psi)^2 ... should always == 1 if psi is normalized
        double V[2*N];          // The potential field
        double _mass;           // The mass of the particle
        double _dt_size;        // The amount of time in nanoseconds that each step represents
        int    _dt_iterations;  // The number of times to apply the time evolution operator
        double _x_min, _x_max;  // The left and right boundaries on the x axis
        double _x_range;        // The distance between the left and right boundaries

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
        void set_dt_size(double dt_size);
        void set_dt_iterations(double dt_iterations);
        void set_x_min(double x_min);
        void set_x_max(double x_max);

        // Getters
    public:
        double const *get_psi() const;
        double const *get_psi_abs2() const;
        double const *get_V() const;
        double get_psi_norm() const;
        double mass() const;
        double dt_size() const;
        double dt_iterations() const;
        double x_min() const;
        double x_max() const;
        double x_range() const;
        double K_max() const;
        double E_min() const;
        double E_max() const;
        double E_range() const;
    };
}

#endif