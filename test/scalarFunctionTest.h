#ifndef SCALAR_FUNCTION_TEST
#define SCALAR_FUNCTION_TEST
/*
 *TESTS FOR SCALAR FUNCTIONS
 */
#include <ScalarValuedFunction.h>
#include <ParametrizedFunction.h>
#include <vector>
#include <MyTypes.h>
#include <cmath>
#include <Curve.h>
namespace FunctionTests{

    static double tolerance = 1E-9;

    /*
     *TEST FOR EXPLICIT FUNCTION
     */
    BOOST_AUTO_TEST_SUITE(Scalar)

    BOOST_AUTO_TEST_SUITE(ExplicitScalarTest)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Scalar::ExplicitScalar::ConstructorTest" << endMsg;
        ExplicitScalarFunction_1D function([](double t) {return std::complex<double>(2*t, 0);});
        msg(1) << "end Scalar::ExplicitScalar::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaulationTest)
    {
        msg(1) << "start Scalar::ExplicitScalar::EvaluationTest" << endMsg;
        ExplicitScalarFunction_1D function([](double t) {return 2*t;});
        BOOST_CHECK_CLOSE(function(1).real(), 2, tolerance);
        BOOST_CHECK_CLOSE(function(3.14).real(), 6.28, tolerance);
        BOOST_CHECK_CLOSE(function(-3.14).real(), -6.28, tolerance);
        msg(1) << "end Scalar::ExplicitScalar::EvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(BoundaryScalarTraceTest)


    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Scalar::BoundaryScalarTraceTest::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        BoundaryScalarTrace_1D function([](double t) {return std::complex<double>(2*t, 0);}, curve);
        msg(1) << "end Scalar::BoundaryScalarTraceTest::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaulationTest)
    {
        msg(1) << "start Scalar::BoundaryScalarTraceTest::EvaluationTest" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        BoundaryScalarTrace_1D function([](double t)->BEM::Complex {return 2*t;}, curve);
        BOOST_CHECK_CLOSE(function(1).real(), 2, tolerance);
        BOOST_CHECK_CLOSE(function(3.14).real(), 6.28, tolerance);
        BOOST_CHECK_CLOSE(function(-3.14).real(), -6.28, tolerance);
        msg(1) << "end Scalar::BoundaryScalarTraceTest::EvaluationTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(TrigonometricParametricFunctionTest)
    
    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Scalar::TrigonometricParametricFunctionTest::ConstructorTest" << endMsg;
        TrigonometricFunction_1D function(1, std::vector<BEM::Complex>{1.0, 0.0});
        msg(1) << "end Scalar::TrigonometricParametricFunctionTest::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ConstantEvaluationTest)
    {
        msg(1) << "start Scalar::TrigonometricParametricFunctionTest::ConstantEvaluationTest" << endMsg;
        TrigonometricFunction_1D function(1, std::vector<BEM::Complex>{1.0, 0.0});
        BOOST_CHECK_CLOSE(std::real(function(0.0)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.1)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.3)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.8)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(1.2)), 1., tolerance);
        msg(1) << "end Scalar::TrigonometricParametricFunctionTest::ConstantEvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest1)
    {
        msg(1) << "start Scalar::TrigonometricParametricFunctionTest::EvaluationTest1" << endMsg;
        TrigonometricFunction_1D function(1., std::vector<BEM::Complex>{.0, 1.});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)),std::sin(2.*M_PI*val), tolerance);
        }
        msg(1) << "end Scalar::TrigonometricParametricFunctionTest::EvaluationTest1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest2)
    {
        msg(1) << "start Scalar::TrigonometricParametricFunctionTest::EvaluationTest2" << endMsg;
        TrigonometricFunction_1D function(1., std::vector<BEM::Complex>{.0, 0., 1.});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)),std::cos(2.*M_PI*val), tolerance);
        }
        msg(1) << "end Scalar::TrigonometricParametricFunctionTest::EvaluationTest2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest3)
    {
        msg(1) << "start Scalar::TrigonometricParametricFunctionTest::EvaluationTest3" << endMsg;
        TrigonometricFunction_1D function(1., std::vector<BEM::Complex>{-.3, 2., 1., 0.5, -0.8});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)),-0.3 + 2*std::sin(2*M_PI*val) + std::cos(2.*M_PI*val) + .5*std::sin(4*M_PI*val) - .8*std::cos(4.*M_PI*val), tolerance);
        }
        msg(1) << "end Scalar::TrigonometricParametricFunctionTest::EvaluationTest3" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(PolynomialParametricFunctionTest)
    
    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Scalar::PolynomialParametricFunctionTest::ConstructorTest" << endMsg;
        PolynomialFunction_1D function(std::vector<BEM::Complex>{1.0, 0.0});
        msg(1) << "end Scalar::PolynomialParametricFunctionTest::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ConstantEvaluationTest)
    {
        msg(1) << "start Scalar::PolynomialParametricFunctionTest::EvaluationTest" << endMsg;
        PolynomialFunction_1D function(std::vector<BEM::Complex>{1.0, 0.0});
        BOOST_CHECK_CLOSE(std::real(function(0.0)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.1)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.3)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.8)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(1.2)), 1., tolerance);
        msg(1) << "end Scalar::PolynomialParametricFunctionTest::EvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest1)
    {
        msg(1) << "start Scalar::PolynomialParametricFunctionTest::EvaluationTest1" << endMsg;
        PolynomialFunction_1D function(std::vector<BEM::Complex>{.0, 1.});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)), val, tolerance);
        }
        msg(1) << "end Scalar::PolynomialParametricFunctionTest::EvaluationTest1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest2)
    {
        msg(1) << "start Scalar::PolynomialParametricFunctionTest::EvaluationTest2" << endMsg;
        PolynomialFunction_1D function(std::vector<BEM::Complex>{.0, 0., 1.});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)),val*val, tolerance);
        }
        msg(1) << "end Scalar::PolynomialParametricFunctionTest::EvaluationTest2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest3)
    {
        msg(1) << "start Scalar::PolynomialParametricFunctionTest::EvaluationTest3" << endMsg;
        PolynomialFunction_1D function(std::vector<BEM::Complex>{-.3, 2., 1., 0.5, -0.8});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)),-0.3 + 2.*val + val*val + 0.5*val*val*val - 0.8*val*val*val*val, tolerance);
        }
        msg(1) << "end Scalar::PolynomialParametricFunctionTest::EvaluationTest3" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()


    BOOST_AUTO_TEST_SUITE(PwConstantFunctionTest)
    
    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Scalar::PwConstantFunctionTest::ConstructorTest" << endMsg;
        PwConstantFunction_1D function(0, 1, std::vector<BEM::Complex>{1.0, 0.0});
        msg(1) << "end Scalar::PwConstantFunctionTest::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ConstantEvaluationTest)
    {
        msg(1) << "start Scalar::PwConstantFunctionTest::EvaluationTest" << endMsg;
        PwConstantFunction_1D function(0, 1, std::vector<BEM::Complex>{1.0, 0.0});
        BOOST_CHECK_CLOSE(std::real(function(0.0)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.1)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.3)), 1., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(.8)), 0., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(1.2)), 0., tolerance);
        msg(1) << "end Scalar::PwConstantFunctionTest::EvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest1)
    {
        msg(1) << "start Scalar::PwConstantFunctionTest::EvaluationTest1" << endMsg;
        PwConstantFunction_1D function(0, 1, std::vector<BEM::Complex>{.0, 1.});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)), (val < 0.5) ? 0.0 : ((val <= 1) ? 1.0 : 0.0), tolerance);
        }
        msg(1) << "end Scalar::PwConstantFunctionTest::EvaluationTest1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest2)
    {
        msg(1) << "start Scalar::PwConstantFunctionTest::EvaluationTest2" << endMsg;
        PwConstantFunction_1D function(0, 1, std::vector<BEM::Complex>{.0, 0., 1.});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        for (const auto &val : values) {
            BOOST_CHECK_CLOSE(std::real(function(val)),(val < 2./3.) ? 0.0 : ((val <= 1.) ? 1.0 : 0.0), tolerance);
        }
        msg(1) << "end Scalar::PwConstantFunctionTest::EvaluationTest2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaluationTest3)
    {
        msg(1) << "start Scalar::PwConstantFunctionTest::EvaluationTest3" << endMsg;
        PwConstantFunction_1D function(0, 1,std::vector<BEM::Complex>{-.3, 2., 1., 0.5, -0.8});
        std::vector<double> values{0.0, 0.1, 0.3, 0.8, 1.2};
        std::vector<double> bp{0.0, 0.2, 0.4, 0.6, .8, 1.};

        BOOST_CHECK_CLOSE(std::real(function(0.0)),-.3, tolerance);
        BOOST_CHECK_CLOSE(std::real(function(0.1)),-.3, tolerance);
        BOOST_CHECK_CLOSE(std::real(function(0.3)),2., tolerance);
        BOOST_CHECK_CLOSE(std::real(function(0.8)),-.8, tolerance);
        BOOST_CHECK_CLOSE(std::real(function(1.2)),0.0, tolerance);
        
        msg(1) << "end Scalar::PwConstantFunctionTest::EvaluationTest3" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()
    
    BOOST_AUTO_TEST_SUITE(BoundaryScalarRestrictionTest)


    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Scalar::BoundaryScalarRestrictionTest::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        BoundaryScalarRestriction_1D function([](double t, double s) ->BEM::Complex {return std::complex<double>(2*t, s);}, curve);
        msg(1) << "end Scalar::BoundaryScalarRestrictionTest::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(EvaulationTest)
    {
        msg(1) << "start Scalar::BoundaryScalarRestrictionTest::EvaluationTest" << endMsg;
        TrigonometricCurve curve(3, 0, {1}, {0});
        BoundaryScalarRestriction_1D function([](const double t, const double s) {return std::complex<double>(2*t, s);}, curve);

        BOOST_CHECK_CLOSE(function(1).real(), 6, tolerance);
        BOOST_CHECK_CLOSE(function(1).imag(), std::sin(2.0*M_PI), tolerance);

        BOOST_CHECK_CLOSE(function(M_PI - 3).real(), 6*M_PI - 18, tolerance);
        BOOST_CHECK_CLOSE(function(M_PI - 3).imag(), std::sin(2.0*M_PI*M_PI), tolerance);
        msg(1) << "end Scalar::BoundaryScalarRestrictionTest::EvaluationTest" << endMsg;
    }


    BOOST_AUTO_TEST_SUITE_END()
    
    BOOST_AUTO_TEST_SUITE_END()

}
#endif
