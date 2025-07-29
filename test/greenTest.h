#ifndef GREEN_TEST
#define GREEN_TEST
#include <cmath>
#include <GreenH.h>
#include <RegularP0.h>
#include <GreenQP.h>
#include <Point.h>
#include <iostream>
#include <WeaklySingular.h>
#include <DiscreteSpace.h>
#include <Curve.h>
#include <TwoDimensionalIntegration.h>
#include <ScalarValuedFunction.h>
#include <EmpiricalInterpolation.h>
#include <Utilities.h>
#include <memory>
#include <chrono>
#include <Mesh.h>
#include <RegularP0Mesh.h>
namespace GreenTest{
    static double tolerance = 1E-6;

    BOOST_AUTO_TEST_SUITE(Green)

    BOOST_AUTO_TEST_SUITE(Helmholtz)

    BOOST_AUTO_TEST_CASE(SingularityIntegration)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegration" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        RegularP0_1D piecewiseConstant(curve, 1);
        GreenH2D green(1);
        green.set1DQuad(40);
        auto &testFunction = piecewiseConstant.basisFunction(0);
        auto &trialFunction = piecewiseConstant.basisFunction(0);
        msg(5) << "Got Base functions" << std::endl;
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -1.5, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, 1).real(), std::log(4.)-1.5, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, -1).real(), std::log(4.)-1.5, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegration" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SingularityIntegrationMesh)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegration" << endMsg;
        GreenH2D green(1);
        TrigonometricCurve curve(1, 0, {0}, {0});
        RegularP0_1D piecewiseConstant(curve, 1);
        StraightCurve straight(1);
        MeshCurve1D mesh(1, straight);
        RegularP0Mesh_1D space(mesh);
        green.set1DQuad(40);
        auto &testFunction = piecewiseConstant.basisFunction(0);
        auto &trialFunction = piecewiseConstant.basisFunction(0);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -1.5, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegration" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(SingularityIntegrationP1)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegrationP1" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        RegularP1_1D piecewiseLinear(curve, 2);
        msg(5) << "Constructed space" << endMsg;
        GreenH2D green(1);
        green.setHighPrecision(1);
        auto &testFunction = piecewiseLinear.basisFunction(1);
        auto &trialFunction = piecewiseLinear.basisFunction(1);
        msg(5) << "Got test and trial function" << endMsg;
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -0.4375, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, 0.0001).real(), -0.4375, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, 1).real(), -0.017068546, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, -1).real(), -0.017068546, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegrationP1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SingularityIntegration2)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegration2" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        RegularP0_1D piecewiseConstant(curve, 10);
        GreenH2D green(1);
        auto &testFunction = piecewiseConstant.basisFunction(0);
        auto &trialFunction = piecewiseConstant.basisFunction(4);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -0.0092156, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, 1.).real(), -0.00513153, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, -1.).real(), 0.00336047, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegration2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SingularityIntegrationP12)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegrationP12" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        RegularP1_1D piecewiseP1(curve, 3);
        RegularP0_1D piecewiseConstant(curve, 1);
        GreenH2D green(1);
        green.setHighPrecision(1);
        auto &testFunction = piecewiseP1.basisFunction(1);
        auto &trialFunction = piecewiseConstant.basisFunction(0);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -0.801142, 1E-3);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, trialFunction, testFunction).real(), -0.801142, 1E-3);
        // BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, 1.).real(), -0.00513153, 1E-3);
        // BOOST_CHECK_CLOSE(2*M_PI*green.integrateShiftedSingularity(curve, testFunction, trialFunction, -1.).real(), 0.00336047, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegrationP12" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(SingularityIntegration3)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegration3" << endMsg;
        TrigonometricCurve curve(1, 0, {0}, {0});
        RegularP0_1D piecewiseConstant(curve, 10);
        RegularP0_1D piecewiseConstant2(curve, 1);
        GreenH2D green(1);
        auto &testFunction = piecewiseConstant.basisFunction(4);
        auto &trialFunction = piecewiseConstant.basisFunction(0);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -0.0092156, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegration3" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(SingularityIntegration4)
    {
        msg(1) << "start Green::Helmholtz::SingularityIntegration4" << endMsg;
        TrigonometricCurve curve(1, 0, {0.8}, {0});
        RegularP0_1D piecewiseConstant(curve, 10);
        RegularP0_1D piecewiseConstant2(curve, 1);
        GreenH2D green(1);
        GaussLegendre_2D<20,20> gauss;
        auto &testFunction = piecewiseConstant.basisFunction(4);
        auto &trialFunction = piecewiseConstant.basisFunction(0);
        BOOST_CHECK_CLOSE(2*M_PI*green.integrateSingularity(curve, testFunction, trialFunction).real(), -0.0092156, 1E-3);
        msg(1) << "end Green::Helmholtz::SingularityIntegration4" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(PECGratingProb1)
    {
        TrigonometricCurve curve(1.0, 0, {0.1}, {-0.2});
        GreenHQP2D greenH(1., 0.15, 2.);
        RegularP0_1D pcSpace(curve, 10);
        DiscreteSpaceOnCurve_1D &space = pcSpace;
        GreenQP2D &green = greenH;
        green.setWindowTerms(30);
        WeaklySingular ws(space, space, green);
        ws.assembleMassMatrix();
        auto pw = BEM::generatePlaneWave(0.15, 2., curve);
        auto rhs = space.testAgainstBasis(*pw);
        BEM::Matrix mat = ws.getMatrix();
        BEM::ColVector sol = mat.colPivHouseholderQr().solve(rhs);
        PECGrating grt(static_cast<PeriodicCurve*>(new TrigonometricCurve(curve)), 0.15, 2., 10, 30);
        Problem &g = grt;
        g.buildDiscrete();
        g.solve();
        BOOST_CHECK_CLOSE((grt.getSolutionVec() - sol).norm(), 0.0, tolerance);
    }
    
    BOOST_AUTO_TEST_CASE(PECGratingProb2)
    {
        TrigonometricCurve curve(1.5, 0, {0.1, 0.02, -0.1}, {-0.2});
        GreenHQP2D greenH(1.5, 0.15, 2.);
        RegularP0_1D pcSpace(curve, 10);
        DiscreteSpaceOnCurve_1D &space = pcSpace;
        GreenQP2D &green = greenH;
        green.setWindowTerms(30);
        WeaklySingular ws(space, space, green);
        ws.assembleMassMatrix();
        auto pw = BEM::generatePlaneWave(0.15, 2., curve);
        auto rhs = space.testAgainstBasis(*pw);
        BEM::Matrix mat = ws.getMatrix();        
        BEM::ColVector sol = mat.colPivHouseholderQr().solve(rhs);
        PECGrating grt(static_cast<PeriodicCurve*>(new TrigonometricCurve(curve)), 0.15, 2., 10, 30);
        Problem &g = grt;
        g.buildDiscrete();
        g.solve();
        BOOST_CHECK_CLOSE((grt.getSolutionVec() - sol).norm(), 0.0, tolerance);
    }


    BOOST_AUTO_TEST_CASE(PECGratingProbMesh1)
    {
        int elements = 200;
        double wavenumber = 6.;
        double angle = -0.15;
        TrigonometricCurve curve(1.0, 0, {-0.1}, {-0.2});
        GreenHQP2D greenH(1., angle, wavenumber);
        GreenQP2D &green = greenH;
        green.setWindowTerms(30);
        PECGrating grt(static_cast<PeriodicCurve*>(new TrigonometricCurve(curve)), angle, wavenumber, elements, 30);
        Problem &g = grt;
        g.buildDiscrete();
        g.solve();

        StraightCurve straight(1.0);
        MeshCurve1D mesh(elements, straight);
        auto transformation = [](Point2D point) {
            double t = point.getX();
            return Point2D(t, -0.1*std::sin(2.0*M_PI*t)-0.2*std::cos(2.0*M_PI*t));
        };
        auto jacobian = [](double t, double s) {
            return Vector2D(-0.1*std::cos(2.0*M_PI*t)+0.2*std::sin(2.0*M_PI*t), -1.).norm();
        };
        auto pw = [wavenumber, angle, &transformation](double t, double s) -> BEM::Complex {
            auto np = transformation(Point2D(t, s));
            double normExp = (wavenumber*(std::sin(angle)*np.getX() - std::cos(angle)*np.getY()));
            return BEM::Complex(std::cos(normExp), std::sin(normExp));
        };
        auto pws = [wavenumber, angle, &transformation](double t) -> BEM::Complex {
            auto np = transformation(Point2D(t, 0));
            double normExp = (wavenumber*(std::sin(angle)*np.getX() - std::cos(angle)*np.getY()));
            return BEM::Complex(std::cos(normExp), std::sin(normExp));
        };
        auto g_pw = BEM::generatePlaneWave(angle, wavenumber, curve);
        BEM::plotFunction("pw", ExplicitScalarFunction_1D(pws));
        BEM::plotFunction("gpw",*g_pw);
        RegularP0Mesh_1D spaceMesh(mesh);
        WeaklySingularMesh ws(spaceMesh, spaceMesh, green, Transformation_2D(std::move(transformation)), ExplicitScalarFunction_2D(std::move(jacobian)));
        ws.assembleMassMatrix();
        BEM::Matrix mat = ws.getMatrix();
        auto rhs = spaceMesh.testAgainstBasis(ExplicitScalarFunction_2D(std::move(pw)));
        BEM::ColVector sol = mat.colPivHouseholderQr().solve(rhs);
        BEM::plotFunction("new", *(g.getSpace().generateFunction(sol)));
        BEM::plotFunction("prob", *(g.getSpace().generateFunction(g.getSolutionVec())));
    }

    
    BOOST_AUTO_TEST_CASE(IntegrationWeaklySingular)
    {
        msg(1) << "start Green::Helmholtz::IntegrationWeaklySingular" << endMsg;
        GreenH2D green(2);
        green.setHighPrecision();
        TrigonometricCurve curve(1.5, 0, {0.1, 0.3}, {0.5});
        RegularP0_1D piecewiseConstant(curve, 200);
        WeaklySingular ws(piecewiseConstant, piecewiseConstant, green);
        ws.setHighPrecision();
        GaussLegendre_2D<200, 201> integ;
        ExplicitScalarFunction_2D g([&](double t, double s) {return green(Point2D(curve.at(t).getX(), curve.at(t).getY()), Point2D(curve.at(s).getX(), curve.at(s).getY()));});
        for (int i = 0; i < 10; ++i) {
            auto &testFunction = piecewiseConstant.basisFunction(i);
            auto testSupport = testFunction.support();
            for (int j = 0; j < 10; ++j) {
                auto &trialFunction = piecewiseConstant.basisFunction(j);
                auto trialSupport = trialFunction.support();
                if (j == i) {
                    continue;
                }
                BOOST_CHECK_CLOSE(ws.indexedDuality(i, j).real(), integ.integrate(testSupport, trialSupport, g).real(), 1E-1);
            }
        }
        msg(1) << "end Green::Helmholtz::IntegrationWeaklySingular" << endMsg;
    }
    
    BOOST_AUTO_TEST_SUITE_END()
    
    BOOST_AUTO_TEST_SUITE(QuasiPer)

    BOOST_AUTO_TEST_CASE(Evaluation)
    {
        msg(1) << "start Green::QuasiPer::Evaluation" << endMsg;
        GreenHQP2D green(1, 0, 1.3);
        Point2D X(0,0), Y(.1,.1);
        green.setWindowTerms(200);
        BOOST_CHECK_CLOSE(green.spectralSum(X, Y).real(), green.windowedSum(X, Y).real(), 1E-4);
        msg(1) << "end Green::QuasiPer::Evaluation" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(QuasiPeriodicity)
    {
        msg(1) << "start Green::QuasiPer::QuasiPeriodicity" << endMsg;
        GreenHQP2D green(1, 5, 1.2);
        auto qp = green.getQP();
        Point2D X(0,0), Y(.1,.1);
        Point2D X1(1,0), Y1(.1,.1);
        BOOST_CHECK_CLOSE(green(X1, Y).real(), (qp*green(X, Y)).real(), 1E-4);
        msg(1) << "end Green::QuasiPer::QuasiPeriodicity" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(WeakMatrix)
    {
        GreenHQP2D green(20, 2, 1, 40);
        GreenHQP2D greenM(20, 2, 1, 40);
        greenM.setHighPrecision();

        TrigonometricCurve curve(1.5, 0, {0.1, 0.3}, {0.5});
        StraightCurve straight(1.);
        MeshCurve1D mesh(20, straight);
        RegularP0Mesh_1D space(mesh);
        RegularP0_1D piecewiseConstant(curve, 20);
        WeaklySingular ws(piecewiseConstant, piecewiseConstant, green);
        auto transformation = [](Point2D point) {
            double t = 1.5*point.getX();
            return Point2D(t, point.getY() + 0.1*std::sin(2.0*M_PI*t/(1.5))+0.3*std::sin(2.0*M_PI*2*t/(1.5))+0.5*std::cos(2.0*M_PI*t/(1.5)));
        };
        WeaklySingularMesh wsm(space, space, greenM, Transformation_2D(std::move(transformation)), ExplicitScalarFunction_2D([](double t, double s){return Point2D(2.*M_PI/(1.5)*(0.1*std::cos(2.0*M_PI*t/(1.5))+0.3*2.*std::cos(2.0*M_PI*2*t/(1.5))-0.5*std::sin(2.0*M_PI*t/(1.5))), 1.5).norm();}));
        wsm.assembleMassMatrix();
        ws.assembleMassMatrix();
        msg(0) << "Rel error: " << (wsm.getMatrix() -  ws.getMatrix()).norm()/ws.getMatrix().norm() << endMsg;
        BOOST_CHECK((wsm.getMatrix() -  ws.getMatrix()).norm()/ws.getMatrix().norm() <  1e-4);
    }
    
    BOOST_AUTO_TEST_CASE(IntegrationWeaklySingular)
    {
        msg(1) << "start Green::QuasiPer::IntegrationWeaklySingular" << endMsg;
        GreenHQP2D green(20, 2, 1, 40);
        TrigonometricCurve curve(1, 0, {0.1}, {0});
        RegularP0_1D piecewiseConstant(curve, 5);
        WeaklySingular ws(piecewiseConstant, piecewiseConstant, green);
        ws.setHighPrecision();
        GaussLegendre_2D<60, 61> integ;
        ExplicitScalarFunction_2D g([&](double t, double s) {return green(Point2D(curve.at(t).getX(), curve.at(t).getY()), Point2D(curve.at(s).getX(), curve.at(s).getY()));});
        for (int i = 0; i < 5; ++i) {
            auto &testFunction = piecewiseConstant.basisFunction(i);
            auto testSupport = testFunction.support();
            for (int j = 0; j < 5; ++j) {
                auto &trialFunction = piecewiseConstant.basisFunction(j);
                auto trialSupport = trialFunction.support();
                if (j == i) {
                    continue;
                }
                BOOST_CHECK_CLOSE(ws.indexedDuality(i, j).real(), integ.integrate(testSupport, trialSupport, g).real(), 1E-2);
            }
        }
        msg(1) << "end Green::QuasiPer::IntegrationWeaklySingular" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(MatrixNormTest, * boost::unit_test::disabled())
    {
        msg(1) << "start Green::QuasiPer::PolynomialIntegrationTest" << endMsg;
        GreenHQP2D green(20, 2, 1, 40);
        TrigonometricCurve curve(1, 0, {0.1}, {0});
        for (int i = 0; i < 3; ++i) {
            RegularP0_1D piecewiseConstant(curve, 10*std::pow(2, i));
            WeaklySingular ws(piecewiseConstant, piecewiseConstant, green);
            ws.assembleMassMatrix();
        }
        msg(1) << "end Green::QuasiPer::PolynomialIntegrationTest" << endMsg;
    }    
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE_END()
    
}

#endif
