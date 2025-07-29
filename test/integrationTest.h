#ifndef ONE_DIM_INTEGRATION_TEST
#define ONE_DIM_INTEGRATION_TEST

#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <ScalarValuedFunction.h>
#include <cmath>
#include <iostream>
#include <boost/math/quadrature/tanh_sinh.hpp>
/*
 *TESTS FOR INTEGRATION
 */

namespace integrationTests{

    static double tolerance = 1E-9;

    BOOST_AUTO_TEST_SUITE(Integration)

    /*
     *TEST FOR 1D tanh_sinh INTEGRATION
     */
    BOOST_AUTO_TEST_SUITE(TanHSinH)
    
    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        boost::math::quadrature::tanh_sinh<double> integrator{};
    }

    BOOST_AUTO_TEST_CASE(IntegrationTest)
    {
        boost::math::quadrature::tanh_sinh<double> integrator{};
        auto res1 = integrator.integrate([](double x) {return 1.0/std::sqrt(x);}, 0, 1);
        auto res2 = integrator.integrate([](double x, double t) {
            return x < 0.25 ? 1.0/std::sqrt(-t) : 1.0/std::sqrt(x);}, 0, 1);
        BOOST_CHECK_CLOSE(res1, res2, tolerance);
    }

    BOOST_AUTO_TEST_SUITE_END()

    
    /*
     *TEST FOR 1D TRAPEZOIDAL INTEGRATION
     */
    BOOST_AUTO_TEST_SUITE(Trapezoidal)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Integration::Trapezoidal::ConstructorTest" << endMsg;
        TrapezoidalQuadrature_1D quad(10);
        BOOST_REQUIRE_THROW(TrapezoidalQuadrature_1D quad2(0), std::invalid_argument);
        BOOST_REQUIRE_THROW(TrapezoidalQuadrature_1D quad3(1), std::invalid_argument);
        msg(1) << "end Integration::Trapezoidal::ConstructorTest" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(ConstantIntegrationTest)
    {
        msg(1) << "start Integration::Trapezoidal::ConstantIntegrationTest" << endMsg;
        TrapezoidalQuadrature_1D quad(2);
        ExplicitScalarFunction_1D function([](double t) {return 1.0 + t*0;});
        BOOST_CHECK_CLOSE(quad.integrate(2, 3, function).real(), 1, tolerance);

        ExplicitScalarFunction_1D function2([](double t) {return 3.0 + t*0;});
        BOOST_CHECK_CLOSE(quad.integrate(-2, 30, function2).real(), 32*3, tolerance);

        ExplicitScalarFunction_1D function3([](double t) {return 1.5 + t*0;});
        BOOST_CHECK_CLOSE(quad.integrate(-2, 0, function3).real(), 3, tolerance);
        msg(1) << "end Integration::Trapezoidal::ConstantIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(AffineIntegrationTest)
    {
        msg(1) << "start Integration::Trapezoidal::AffineIntegrationTest" << endMsg;
        TrapezoidalQuadrature_1D quad(2);
        ExplicitScalarFunction_1D function([](double t) {return 1.0 + t;});
        BOOST_CHECK_CLOSE(quad.integrate(2, 3, function).real(), 3.5, tolerance);
        BOOST_CHECK_CLOSE(quad.integrate(-1, 3, function).real(), 8, tolerance);
        msg(1) << "end Integration::Trapezoidal::AffineIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TrigonometricIntegrationTest)
    {
        msg(1) << "start Integration::Trapezoidal::TrigonometricIntegrationTest" << endMsg;
        TrapezoidalQuadrature_1D quad(1000);
        ExplicitScalarFunction_1D function([](double t) {return std::sin(2*M_PI*t) + 0.1;});
        BOOST_CHECK_CLOSE(quad.integrate(-1, 1, function).real(), 0.2, tolerance);
        BOOST_CHECK_CLOSE(quad.integrate(-1, 3, function).real(), 0.4, tolerance);
        
        ExplicitScalarFunction_1D function2([](double t) {return std::cos(2*M_PI*t);});
        BOOST_CHECK_CLOSE(1.0/M_PI, quad.integrate(-0.25, 0.25, function2).real(), 1E-4);
        msg(1) << "end Integration::Trapezoidal::TrigonometricIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()

    /*
     *TEST FOR 1D GAUSS-LEGENDRE INTEGRATION
     */
    BOOST_AUTO_TEST_SUITE(GaussLegendre)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Integration::GaussLegendre::ConstructorTest" << endMsg;
        GaussLegendre_1D quad(10);
        msg(1) << "end Integration::GaussLegendre::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ConstantIntegrationTest)
    {
        msg(1) << "start Integration::GaussLegendre::ConstantIntegrationTeset" << endMsg;
        GaussLegendre_1D quad(7);
        ExplicitScalarFunction_1D function([](double t) {return 1.0 + t*0;});
        BOOST_CHECK_CLOSE(quad.integrate(2, 3, function).real(), 1, tolerance);

        ExplicitScalarFunction_1D function2([](double t) {return 3.0 + t*0;});
        BOOST_CHECK_CLOSE(quad.integrate(-2, 30, function2).real(), 32*3, tolerance);

        ExplicitScalarFunction_1D function3([](double t) {return 1.5 + t*0;});
        BOOST_CHECK_CLOSE(quad.integrate(-2, 0, function3).real(), 3, tolerance);
        msg(1) << "end Integration::GaussLegendre::ConstantIntegrationTeset" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(AffineIntegrationTest)
    {
        msg(1) << "start Integration::GaussLegendre::AffineIntegrationTest" << endMsg;
        GaussLegendre_1D quad(3);
        ExplicitScalarFunction_1D function([](double t) {return 1.0 + t;});
        BOOST_CHECK_CLOSE(quad.integrate(2, 3, function).real(), 3.5, tolerance);
        BOOST_CHECK_CLOSE(quad.integrate(-1, 3, function).real(), 8, tolerance);
        msg(1) << "end Integration::GaussLegendre::AffineIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TrigonometricIntegrationTest)
    {
        msg(1) << "start Integration::GaussLegendre::TrigonometricIntegrationTest" << endMsg;
        GaussLegendre_1D quad(13);
        ExplicitScalarFunction_1D function([](double t) {return std::sin(2*M_PI*t) + 0.1;});
        BOOST_CHECK_CLOSE(quad.integrate(-1, 1, function).real(), 0.2, tolerance);
        BOOST_CHECK_CLOSE(quad.integrate(-1, 3, function).real(), 0.4, tolerance);
        
        ExplicitScalarFunction_1D function2([](double t) {return std::cos(2*M_PI*t);});
        BOOST_CHECK_CLOSE(1.0/M_PI, quad.integrate(-0.25, 0.25, function2).real(), 1E-4);
        msg(1) << "end Integration::GaussLegendre::TrigonometricIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(GaussLegendre2D)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Integration::GaussLegendre2D::ConstructorTest" << endMsg;
        GaussLegendre_2D<10, 12> quad;
        msg(1) << "end Integration::GaussLegendre2D::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ConstantIntegrationTest)
    {
        msg(1) << "start Integration::GaussLegendre2D::ConstantIntegrationTest" << endMsg;
        GaussLegendre_2D<7, 5> quad;
        ExplicitScalarFunction_2D function([](double t, double s) {return 0.*(t+s) + 1.0;});
        BOOST_CHECK_CLOSE(quad.integrate(2, 3, -1, 3, function).real(), 4, tolerance);

        ExplicitScalarFunction_2D function2([](double t, double s) {return 0.*(t+s) + M_PI;});
        BOOST_CHECK_CLOSE(quad.integrate(0, 0, -1, 3, function2).real(), 0, tolerance);;
        BOOST_CHECK_CLOSE(quad.integrate(-1, 10, -1, 3, function2).real(), M_PI*44, tolerance);

        ExplicitScalarFunction_2D function3([](double t, double s) {return 0.*(t+s) + -2;});
        BOOST_CHECK_CLOSE(quad.integrate(-1, 1, -1, 1, function3).real(), -8, tolerance);
        msg(1) << "end Integration::GaussLegendre2D::ConstantIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(PolynomialIntegrationTest)
    {
        msg(1) << "start Integration::GaussLegendre2D::PolynomialIntegrationTest" << endMsg;
        GaussLegendre_2D<3, 3> quad;
        ExplicitScalarFunction_2D function([](double t, double s) {return 1.0 + 3*t - 2*s;});
        BOOST_CHECK_CLOSE(quad.integrate(2, 3, -1, 3, function).real(), 26, tolerance);

        ExplicitScalarFunction_2D function2([](double t, double s) {return M_PI - 2*t + 3*t*s - s*s;});
        BOOST_CHECK_CLOSE(quad.integrate(-1, 10, -1, 3, function2).real(), M_PI*44 + 286.0/3.0, tolerance);
        msg(1) << "end Integration::GaussLegendre2D::PolynomialIntegrationTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(AdaptiveTrapezoidalTest)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Integration::AdaptiveTrapezoidalTest::ConstructorTest" << endMsg;
        AdaptiveTrapezoidal_1D integrator(1e-2);
        msg(1) << "end Integration::AdaptiveTrapezoidalTest::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(PolynomialIntegrationTest)
    {
        msg(1) << "start Integration::AdaptiveTrapezoidalTest::PolynomialIntegrationTest" << endMsg;
        AdaptiveTrapezoidal_1D integrator(1e-5);
        ExplicitScalarFunction_1D function([](double t) {return 1.0 + 3*t*t;});
        BOOST_CHECK_CLOSE(std::real(integrator.integrate(0, 1, function)), 2, 1e-4);
        msg(1) << "end Integration::AdaptiveTrapezoidalTest::PolynomialIntegrationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SingularIntegrationTest)
    {
        msg(1) << "start Integration::AdaptiveTrapezoidalTest::SingularIntegrationTest" << endMsg;
        AdaptiveTrapezoidal_1D integrator(1e-6);
        ExplicitScalarFunction_1D function([](double t) {return std::pow(t, 0.4);});
        BOOST_CHECK_CLOSE(std::real(integrator.integrate(0, 1, function)), 5./7., 1e-4);
        msg(1) << "end Integration::AdaptiveTrapezoidalTest::SingularIntegrationTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()
    
    BOOST_AUTO_TEST_SUITE_END()
    
    
}

#endif
