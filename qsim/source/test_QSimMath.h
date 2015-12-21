#include "gtest/gtest.h"
#include "qsim/QSimMath.h"
#include "qsim/QSimModel.h"
#include "qsim/WaveFunctionPresets.h"

TEST(MathematicaComparisonTest, TimeEvolution) {

    // The simulation was performed with high precision
    // using Mathematica with specific values. Samples
    // from the results of the time evolution operation
    // have been chosen for comparison.

    // Set up wave functions
    double x_min        = 0.0;
    double x_max        = 100.0;
    double x0           = (x_min + x_max) * 0.3;
    double k0           = 10.0;
    double alpha        = 1.0;
    double V_max        = 8.0;
    double thickness    = 0.48 * (x_max - x_min);
    qsim::Gaussian psi_0(x0, k0, alpha);
    qsim::SquareWell V_0(x_min + thickness, x_max - thickness, V_max);
    qsim::QSimModel model(psi_0, V_0);

    // Set up model
    model.set_dt_size(1.0);
    model.set_dt_iterations(1);

    // Begin making wild assertions!!!
    ASSERT_EQ(qsim::QSimMath::M(&model), 34);
}