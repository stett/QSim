#ifndef QSIMMODEL_H
#define QSIMMODEL_H
#include <gsl/gsl_complex.h>
#include "qsim/QSimConstants.h"
#include "qsim/QSimCoordinates.h"

namespace qsim {
    class QSimModel {

        // Make friends with the QSim math class
        friend class QSimMath;

        // Public types
    public:
        //typedef gsl_complex (*WaveFunction)(double x, double x_min, double x_max);

        // Private members
    private:
        //WaveFunction psi_0;     // The analytical initial wavefunction, a pointer to a function of the form z=f(x)
        //WaveFunction V_0;       // The analytical potential function, a pointer to a function of the form V=f(x)
        double psi[2*N];        // The wave function of a the particle
        double psi_abs2[N];     // Since abs(psi)^2 is used in multiple places, we will just always calculate it once
        double psi_norm;        // The integral of abs(psi)^2 ... should always == 1 if psi is normalized
        double V[2*N];          // The potential field
        double _mass;           // The mass of the particle
        double _dt_size;        // The amount of time in nanoseconds that each step represents
        int    _dt_iterations;  // The number of times to apply the time evolution operator
        double _x_min, _x_max;  // The left and right boundaries on the x axis

        // 'Tors
    public:
        QSimModel();
        ~QSimModel();

        /*
        // Copying
        QSimModel(const QSimModel &model);
        QSimModel &operator=(const QSimModel &model);
        */

        // Public methods
    public:
        void evolve();

        // Internal methods
    private:
        void compute_psi_abs2();

        // Setters
    public:
        template<typename Psi0, typename V0>
        void set_functions(const Psi0 &psi_0, const V0 &V_0);
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

    template<typename Psi0, typename V0>
    void QSimModel::set_functions(const Psi0 &psi_0, const V0 &V_0) {

        // Temp vars
        gsl_complex psi_val;
        gsl_complex V_val;
        double x;

        // Loop through the points of the psi and V arrays & save the calculated values
        for (int n = 0; n < N; ++n) {

            // Find this x-coordinate
            x = QSIM_COORD_INDEX_TO_SPACE_X(n, x_min(), x_range());

            // Get the values of the psi and V functions at this point
            psi_val = psi_0(x);
            V_val   = V_0(x);

            // Record the values
            GSL_COMPLEX_PACKED_REAL(psi, 1, n) = GSL_REAL(psi_val);
            GSL_COMPLEX_PACKED_IMAG(psi, 1, n) = GSL_IMAG(psi_val);
            GSL_COMPLEX_PACKED_REAL(V, 1, n)   = GSL_REAL(V_val);
            GSL_COMPLEX_PACKED_IMAG(V, 1, n)   = GSL_IMAG(V_val);
        }

        // Get the initial psi_abs2
        compute_psi_abs2();
    }
}


#endif