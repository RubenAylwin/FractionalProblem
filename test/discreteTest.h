#ifndef DISCRETE_TEST
#define DISCRETE_TEST
/*
 *TESTS FOR DISCRETE SPACES
 */
#include <PECGrating.h>
#include <DiscreteSpace.h>
#include <Utilities.h>
#include <Curve.h>
#include <vector>
#include <complex>
#include <boost/test/execution_monitor.hpp>
#include <cmath>
#include <stdexcept>
#include <Utilities.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <DiscreteSpaceMatrixMgr.h>
#include <RegularP0.h>
#include <Mesh.h>
#include <RegularP0Mesh.h>
namespace discreteTests{

    static double tolerance = 1E-9;

    BOOST_AUTO_TEST_SUITE(Discrete)
    /*
     *TEST FOR P0 SPACE
     */
    BOOST_AUTO_TEST_SUITE(P0)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P0::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP0_1D space(curve, 10);
        msg(1) << "end Discrete::P0::ConstructorTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P0::EvaluationTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP0_1D space(curve, 10);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(1/curve.jacobian(0.05), firstFun(0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun(0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(9);
        BOOST_CHECK_CLOSE(0, lastFun(0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1/curve.jacobian(0.95), lastFun(0.95).real(), tolerance);
        msg(1) << "end Discrete::P0::EvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P0::InvalidFunctionTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP0_1D space(curve, 10);

        BOOST_REQUIRE_THROW(space.basisFunction(10), std::invalid_argument);
        msg(1) << "end Discrete::P0::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P0::TestingConstantFunction" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP0_1D space(curve, 20);
        ExplicitScalarFunction_1D function([&](double t){ return 1.*curve.jacobian(t);});
        auto testingResult = space.testAgainstBasis(function);
        for (const auto val : testingResult) {
            BOOST_CHECK_CLOSE(val.real(), 1.0/20, tolerance);
        }
        msg(1) << "end Discrete::P0::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingP0DOFs)
    {
        msg(1) << "start Discrete::P0::TestignP0DOFs" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        unsigned sizeOfSpace = 30;
        RegularP0_1D space(curve, sizeOfSpace);
        for (unsigned i = 0; i < sizeOfSpace; ++i) {
            auto &function = space.basisFunction(i);
            for (unsigned j = 0; j < sizeOfSpace; ++j) {
                auto &DoF = space.degreeOfFreedom(j);
                BOOST_CHECK_CLOSE(std::real(DoF(function)), (i == j ? 1 : 0), tolerance);
            }
        }
        msg(1) << "end Discrete::P0::TestignP0DOFs" << endMsg;
    }
    
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(P1)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P1::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 10);
        msg(1) << "end Discrete::P1::ConstructorTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(PlotTest, * boost::unit_test::disabled())
    {
        msg(1) << "start Discrete::P1::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        RegularP1_1D space(curve, 3);
        auto &firstFun = space.basisFunction(1);
        auto antiDerFirstFun = firstFun.antiDerivative(BasisFunction_1D::Direction::RIGHT);
        auto &antiDerFirstFunRef = *antiDerFirstFun;
        ExplicitScalarFunction_1D anti([&](double t)->BEM::Complex {return antiDerFirstFunRef(t,-0.4);});

        msg(1) << "end Discrete::P1::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SupportTest)
    {
        msg(1) << "start Discrete::P1::SupportTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 11);
        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(0, firstFun.support().first, tolerance);
        BOOST_CHECK_CLOSE(0.1, firstFun.support().second, tolerance);

        auto &randFun = space.basisFunction(4);
        BOOST_CHECK_CLOSE(0.3, randFun.support().first, tolerance);
        BOOST_CHECK_CLOSE(0.5, randFun.support().second, tolerance);

        auto &lastFun = space.basisFunction(10);
        BOOST_CHECK_CLOSE(0.9, lastFun.support().first, tolerance);
        BOOST_CHECK_CLOSE(1, lastFun.support().second, tolerance);
        msg(1) << "end Discrete::P1::SupportTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P1::EvaluationTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 11);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(0.5/curve.jacobian(0.05), firstFun(0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun(0.15).real(), tolerance);

        
        auto &lastFun = space.basisFunction(10);
        BOOST_CHECK_CLOSE(0, lastFun(0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.5/curve.jacobian(0.95), lastFun(0.95).real(), tolerance);

        msg(1) << "end Discrete::P1::EvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P1::InvalidFunctionTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 10);

        BOOST_REQUIRE_THROW(space.basisFunction(10), std::invalid_argument);
        msg(1) << "end Discrete::P1::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P1::TestingConstantFunction" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 21);
        ExplicitScalarFunction_1D function([&](double t){ return 1.*curve.jacobian(t);});
        auto testingResult = space.testAgainstBasis(function);
        BOOST_CHECK_CLOSE(testingResult[0].real(), 1.0/40, tolerance);
        BOOST_CHECK_CLOSE(testingResult[20].real(), 1.0/40, tolerance);
        for (int i = 1; i < 20; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), 1.0/20, tolerance);
        }
        msg(1) << "end Discrete::P1::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingP1DOFs)
    {
        msg(1) << "start Discrete::P1::TestignP1DOFs" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        unsigned sizeOfSpace = 31;
        RegularP1_1D space(curve, sizeOfSpace);
        for (unsigned i = 0; i < sizeOfSpace; ++i) {
            auto &function = space.basisFunction(i);
            for (unsigned j = 0; j < sizeOfSpace; ++j) {
                auto &DoF = space.degreeOfFreedom(j);
                BOOST_CHECK_CLOSE(std::real(DoF(function)), (i == j ? 1 : 0), tolerance);
            }
        }
        msg(1) << "end Discrete::P1::TestignP1DOFs" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingP1Derivative)
    {
        msg(1) << "start Discrete::P1::TestingP1Derivative" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        RegularP1_1D space(curve, 9);

        {
            auto &function = space.basisFunction(0);
            auto &derivative = *(function.derivative());
            BOOST_CHECK_CLOSE(std::real(derivative(0.05)), -8, tolerance);
        }

        {
            auto &function = space.basisFunction(8);
            auto &derivative = *(function.derivative());
            BOOST_CHECK_CLOSE(std::real(derivative(0.95)), 8, tolerance);
        }

        for (int i = 1; i < 8; ++i) {
            auto &function = space.basisFunction(i);
            auto &derivative = *(function.derivative());
            auto support = function.brokenSupport();
            int j = 1;
            for (const auto interval : support) {
                auto a = interval.first;
                auto b = interval.second;
                double mult = j == 1 ? 1. : -1.;
                j *= -1;
                double numDer = mult*std::real(1)/(b-a);
                BOOST_CHECK_CLOSE(std::real(derivative(0.5*a + 0.5*b)), numDer, tolerance);
            }
        }
        msg(1) << "end Discrete::P1::TestingP1Derivative" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE(P1_0)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P1_0::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 10, true);
        msg(1) << "end Discrete::P1_0::ConstructorTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(PlotTest, * boost::unit_test::disabled())
    {
        msg(1) << "start Discrete::P1_0::ConstructorTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        RegularP1_1D space(curve, 3, true);
        auto &firstFun = space.basisFunction(2);
        auto antiDerFirstFun = firstFun.antiDerivative(BasisFunction_1D::Direction::RIGHT);
        auto &antiDerFirstFunRef = *antiDerFirstFun;
        ExplicitScalarFunction_1D anti([&](double t)->BEM::Complex {return antiDerFirstFunRef(t,-0.4);});

        msg(1) << "end Discrete::P1_0::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SupportTest)
    {
        msg(1) << "start Discrete::P1_0::SupportTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 9, true);
        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(0, firstFun.support().first, tolerance);
        BOOST_CHECK_CLOSE(0.2, firstFun.support().second, tolerance);

        auto &randFun = space.basisFunction(4);
        BOOST_CHECK_CLOSE(0.4, randFun.support().first, tolerance);
        BOOST_CHECK_CLOSE(0.6, randFun.support().second, tolerance);

        auto &lastFun = space.basisFunction(8);
        BOOST_CHECK_CLOSE(0.8, lastFun.support().first, tolerance);
        BOOST_CHECK_CLOSE(1, lastFun.support().second, tolerance);
        msg(1) << "end Discrete::P1_0::SupportTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P1_0::EvaluationTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 9, true);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(0, firstFun(0.0).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun(0.0).imag(), tolerance);

        BOOST_CHECK_CLOSE(1, firstFun(0.1).real()*curve.jacobian(0.1), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun(0.0).imag(), tolerance);

        BOOST_CHECK_CLOSE(0.5, firstFun(0.05).real()*curve.jacobian(0.05), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun(0.05).imag(), tolerance);

        auto &lastFun = space.basisFunction(8);
        BOOST_CHECK_CLOSE(0, lastFun(1.).real(), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun(1.).imag(), tolerance);

        BOOST_CHECK_CLOSE(1, lastFun(.9).real()*curve.jacobian(.9), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun(0.9).imag(), tolerance);

        BOOST_CHECK_CLOSE(0.5, lastFun(0.95).real()*curve.jacobian(0.95), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun(0.95).imag(), tolerance);


        msg(1) << "end Discrete::P1_0::EvaluationTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P1_0::InvalidFunctionTest" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 10);

        BOOST_REQUIRE_THROW(space.basisFunction(10), std::invalid_argument);
        msg(1) << "end Discrete::P1_0::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P1_0::TestingConstantFunction" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        RegularP1_1D space(curve, 19, true);
        ExplicitScalarFunction_1D function([&](double t){ return 1.*curve.jacobian(t);});
        auto testingResult = space.testAgainstBasis(function);
        for (unsigned i = 0; i < 19; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), 1.0/20, tolerance);
        }
        msg(1) << "end Discrete::P1_0::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingP1_0DOFs)
    {
        msg(1) << "start Discrete::P1_0::TestingP1_0DOFs" << endMsg;
        TrigonometricCurve curve(1, 0, std::vector<double>{1}, std::vector<double>{0});
        unsigned sizeOfSpace = 29;
        RegularP1_1D space(curve, sizeOfSpace, true);
        for (unsigned i = 0; i < sizeOfSpace; ++i) {
            auto &function = space.basisFunction(i);
            double midPoint = 1./30.*(i + 1.);
            BOOST_CHECK_CLOSE(1.0, std::real(function(midPoint))*curve.jacobian(midPoint), tolerance);
            for (unsigned j = 0; j < sizeOfSpace; ++j) {
                auto &DoF = space.degreeOfFreedom(j);
                double midPoint = 1./30.*(j + 1.);
                if (i==j) {
                    BOOST_CHECK_CLOSE(std::real(DoF(function)), 1., tolerance);
                } else {
                    BOOST_CHECK_CLOSE(std::real(DoF(function)), std::real(function(midPoint)), tolerance);
                }
            }
        }
        msg(1) << "end Discrete::P1_0::TestignP1_0DOFs" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE_END()
}
#endif
