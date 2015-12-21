#include <limits>
#include <gtest/gtest.h>
#include <gsl/gsl_complex.h>
#include "qsim/QSimMath.h"
#include "qsim/QSimModel.h"
#include "qsim/WaveFunctionPresets.h"


class MathematicaComparisonTests : public ::testing::Test {
public:

    // The simulation was performed with high precision
    // using Mathematica with specific values. Samples
    // from the results of the time evolution operation
    // have been chosen for comparison.

    double x_min;
    double x_max;
    qsim::QSimModel model;


    MathematicaComparisonTests() : ::testing::Test(), x_min(0.0), x_max(100.0) {

        double x0           = (x_min + x_max) * 0.3;
        double k0           = 10.0;
        double alpha        = 1.5;
        double V_max        = 3.0;
        double thickness    = 0.48 * (x_max - x_min);

        qsim::Gaussian psi_0(x0, k0, alpha);
        qsim::SquareWell V_0(x_min + thickness, x_max - thickness, V_max);

        // Set up model
        model.set_dt_size(1.0);
        model.set_dt_iterations(1);
        model.set_functions(psi_0, V_0);
    }

protected:
    virtual void SetUp() {}
    virtual void TearDown() {}
};

const double min_double = 0.00000000000001;
#define EXPECT_CLOSE(expected, actual) \
    EXPECT_GE(min_double, abs((double)expected - (double)actual))
    //EXPECT_TRUE(expected > actual ? expected - actual <= min_double : actual - expected <= min_double)
    //EXPECT_GE(min_double, abs((double)expected - (double)actual))

TEST_F(MathematicaComparisonTests, TestKMax) {
    EXPECT_CLOSE(9.85740413756665, model.K_max());
}

TEST_F(MathematicaComparisonTests, TestEMin) {
    EXPECT_CLOSE(0.0, model.E_min());
}

TEST_F(MathematicaComparisonTests, TestEMax) {
    EXPECT_CLOSE(12.85740413756665, model.E_max());
}

TEST_F(MathematicaComparisonTests, TestAlpha) {
    EXPECT_CLOSE(9.76691864317238, qsim::QSimMath::alpha(&model));
}

TEST_F(MathematicaComparisonTests, TestLambda) {
    gsl_complex lambda = qsim::QSimMath::lambda(&model);
    EXPECT_CLOSE(0.999469353473798, GSL_REAL(lambda));
    EXPECT_CLOSE(-0.03257317096426748, GSL_IMAG(lambda));
}

TEST_F(MathematicaComparisonTests, TestChebyshevTermCount) {
    EXPECT_EQ(34, qsim::QSimMath::M(&model));
}

TEST_F(MathematicaComparisonTests, TestTimeEvolution) {
    model.evolve();
    EXPECT_CLOSE(0.103977507130599, model.get_psi_abs2()[150]);
    EXPECT_CLOSE(0.999999999999969, model.get_psi_norm());
}