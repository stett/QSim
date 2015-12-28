#include <limits>
#include <gtest/gtest.h>
#include <gsl/gsl_complex.h>
#include "qsim/QSimMath.h"
#include "qsim/QSimModel.h"
#include "qsim/QSimConstants.h"
#include "qsim/QSimCoordinates.h"
#include "qsim/WaveFunctionPresets.h"


class MathematicaComparisonTests : public ::testing::Test {
public:

    // The simulation was performed with high precision
    // using Mathematica with specific values. Samples
    // from the results of the time evolution operation
    // have been chosen for comparison.

    double x_min;
    double x_max;
    double psi[2*qsim::N];
    double phi_0[2*qsim::N];
    double phi_1[2*qsim::N];
    double phi_m[2*qsim::N];
    qsim::Gaussian psi_0;
    qsim::SquareBarrier V_0;
    qsim::QSimModel model;
    gsl_complex val;

    MathematicaComparisonTests() : ::testing::Test(), x_min(0.0), x_max(100.0) {

        double k0           = 10.0;
        double x0 = (x_min + x_max) * 0.3;
        double alpha        = 1.5;
        double V_max        = 3.0;
        double thickness    = 0.48 * (x_max - x_min);

        psi_0 = qsim::Gaussian(x0, k0, alpha);
        V_0 = qsim::SquareBarrier(x_min + thickness, x_max - thickness, V_max);

        // Set up model
        model.set_dt_size(1.0);
        model.set_dt_iterations(1);
        model.set_functions(psi_0, V_0);

        // Copy the initial data
        for (int i = 0; i < 2 * qsim::N; ++i)
            psi[i] = model.get_psi()[i];
    }

protected:
    virtual void SetUp() {}
    virtual void TearDown() {}
};

const double min_double = 0.0000000000001;
#define EXPECT_CLOSE(expected, actual) \
    EXPECT_GE(min_double, abs((double)(expected) - (double)(actual)))
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

TEST_F(MathematicaComparisonTests, TestCoordinateConversions) {
    double space_x = QSIM_COORD_INDEX_TO_SPACE_X(149, model.x_min(), model.x_range());
    EXPECT_CLOSE(29.1015625, space_x);
}

TEST_F(MathematicaComparisonTests, TestPsiFunctionSample) {
    double x = 29.1015625;
    val = psi_0(x);
    EXPECT_CLOSE(-0.4264866663801281, GSL_REAL(val));
    EXPECT_CLOSE(-0.2009916044119778, GSL_IMAG(val));
}

TEST_F(MathematicaComparisonTests, TestPsiDiscreteSample) {
    val = GSL_COMPLEX_PACKED_GET(psi, 1, 149);
    EXPECT_CLOSE(-0.4264866663801281, GSL_REAL(val));
    EXPECT_CLOSE(-0.2009916044119778, GSL_IMAG(val));
}

TEST_F(MathematicaComparisonTests, TestFFTSecondDerivative) {

    // Perform a second derivative on the model's waveform
    qsim::QSimMath::D2(psi, model.x_range());

    // Compare a point to Mathematica results
    val = GSL_COMPLEX_PACKED_GET(psi, 1, 149);
    EXPECT_CLOSE(43.5290118115666, GSL_REAL(val));
    EXPECT_CLOSE(18.43282851671061, GSL_IMAG(val));
}

TEST_F(MathematicaComparisonTests, TestHNorm) {

    // Use the normalized hamiltonian on the model
    qsim::QSimMath::HNorm(psi, &model);

    // Pull out a specific value and test it against the Mathematica results
    val = GSL_COMPLEX_PACKED_GET(psi, 1, 149);
    EXPECT_CLOSE(0.1685112015740519, GSL_REAL(val));
    EXPECT_CLOSE(0.0917491174730588, GSL_IMAG(val));
}

TEST_F(MathematicaComparisonTests, TestChebyshevTermCount) {
    EXPECT_EQ(34, qsim::QSimMath::M(&model));
}

TEST_F(MathematicaComparisonTests, TestChebyshevTermCoefficient) {
    EXPECT_CLOSE(-0.2290728865881845,   qsim::QSimMath::a_m(&model, 0));
    EXPECT_CLOSE(0.2015985052575604,    qsim::QSimMath::a_m(&model, 1));
    EXPECT_CLOSE(0.499427678477892,     qsim::QSimMath::a_m(&model, 2));
    EXPECT_CLOSE(0.002939976826354512,  qsim::QSimMath::a_m(&model, 3));
    EXPECT_CLOSE(0.00001522681397088471,qsim::QSimMath::a_m(&model, 20));
}

TEST_F(MathematicaComparisonTests, TestChebyshevTerms) {

    // Get the initial terms, 0 and 1
    qsim::QSimMath::phi_01(&model, phi_0, phi_1);
    val = GSL_COMPLEX_PACKED_GET(phi_1, 1, 149);
    EXPECT_CLOSE(0.0917491174730588, GSL_REAL(val));
    EXPECT_CLOSE(-0.1685112015740519, GSL_IMAG(val));

    // Test term 2
    qsim::QSimMath::phi_m(&model, phi_m, phi_0, phi_1);
    qsim::QSimMath::phi_01(&model, phi_m, phi_0, phi_1);
    val = GSL_COMPLEX_PACKED_GET(phi_m, 1, 149);
    EXPECT_CLOSE(-0.2911928251639383, GSL_REAL(val));
    EXPECT_CLOSE(-0.117317049914739, GSL_IMAG(val));

    // Test term 3
    qsim::QSimMath::phi_m(&model, phi_m, phi_0, phi_1);
    qsim::QSimMath::phi_01(&model, phi_m, phi_0, phi_1);
    val = GSL_COMPLEX_PACKED_GET(phi_m, 1, 149);
    EXPECT_CLOSE(0.1987484139031981, GSL_REAL(val));
    EXPECT_CLOSE(-0.3953024599711036, GSL_IMAG(val));

    // Test term 33
    for (int ii = 0; ii < 30; ++ii) {
        qsim::QSimMath::phi_m(&model, phi_m, phi_0, phi_1);
        qsim::QSimMath::phi_01(&model, phi_m, phi_0, phi_1);
    }
    val = GSL_COMPLEX_PACKED_GET(phi_m, 1, 149);
    EXPECT_CLOSE(0.04987795692833926, GSL_REAL(val));
    EXPECT_CLOSE(-0.06800783050485359, GSL_IMAG(val));
}

TEST_F(MathematicaComparisonTests, TestComplexArrayMultiplication) {
    qsim::QSimMath::phi_01(&model, phi_0, phi_1);
    qsim::QSimMath::phi_m(&model, phi_m, phi_0, phi_1);
    double a_m = qsim::QSimMath::a_m(&model, 2);
    qsim::QSimMath::multiply_real(phi_m, a_m);
    val = GSL_COMPLEX_PACKED_GET(phi_m, 1, 149);
    EXPECT_CLOSE(-0.1454297566610444, GSL_REAL(val));
    EXPECT_CLOSE(-0.05859138188479309, GSL_IMAG(val));
}

TEST_F(MathematicaComparisonTests, TestComplexArrayAddition) {
    for (int i = 0; i < 2 * qsim::N; ++i) psi[i] = 0;
    qsim::QSimMath::phi_01(&model, phi_0, phi_1);
    qsim::QSimMath::phi_m(&model, phi_m, phi_0, phi_1);
    double a_0 = qsim::QSimMath::a_m(&model, 0);
    double a_1 = qsim::QSimMath::a_m(&model, 1);
    double a_2 = qsim::QSimMath::a_m(&model, 2);

    qsim::QSimMath::multiply_real(phi_0, a_0);
    qsim::QSimMath::add(psi, phi_0);
    val = GSL_COMPLEX_PACKED_GET(psi, 1, 149);
    EXPECT_CLOSE(0.097696531759068, GSL_REAL(val));
    EXPECT_CLOSE(0.04604172700264225, GSL_IMAG(val));

    qsim::QSimMath::multiply_real(phi_1, a_1);
    qsim::QSimMath::add(psi, phi_1);
    val = GSL_COMPLEX_PACKED_GET(psi, 1, 149);
    EXPECT_CLOSE(0.116193016700337, GSL_REAL(val));
    EXPECT_CLOSE(0.01207012064615791, GSL_IMAG(val));

    qsim::QSimMath::multiply_real(phi_m, a_2);
    qsim::QSimMath::add(psi, phi_m);
    val = GSL_COMPLEX_PACKED_GET(psi, 1, 149);
    EXPECT_CLOSE(-0.02923673996070744, GSL_REAL(val));
    EXPECT_CLOSE(-0.04652126123863519, GSL_IMAG(val));
}


TEST_F(MathematicaComparisonTests, TestComplexMultiplication) {

    gsl_complex a, b, ab;
    GSL_REAL(a) = 1.234567890;
    GSL_IMAG(a) = 0.987654321;
    GSL_REAL(b) = 2.345678901;
    GSL_IMAG(b) = 1.098765432;
    ab = gsl_complex_mul(a, b);
    EXPECT_CLOSE(1.8106994247448571, GSL_REAL(ab));
    EXPECT_CLOSE(3.673220423240359, GSL_IMAG(ab));

    gsl_complex lambda = qsim::QSimMath::lambda(&model);
    gsl_complex al = gsl_complex_mul(a, lambda);
    EXPECT_CLOSE(1.2660838038893418, GSL_REAL(al));
    EXPECT_CLOSE(0.9469164347175083, GSL_IMAG(al));
}


TEST_F(MathematicaComparisonTests, TestTimeEvolution) {
    qsim::QSimMath::U(&model);
    val = GSL_COMPLEX_PACKED_GET(model.get_psi(), 1, 149);

    // With final multiplication by lambda
    EXPECT_CLOSE(0.1035163405993156, GSL_REAL(val));
    EXPECT_CLOSE(0.3053880717374624, GSL_IMAG(val));
}


TEST_F(MathematicaComparisonTests, TestNormalization) {
    model.evolve();
    EXPECT_CLOSE(0.103977507130599, model.get_psi_abs2()[149]);
    EXPECT_CLOSE(0.999999999999969, model.get_psi_norm());
}