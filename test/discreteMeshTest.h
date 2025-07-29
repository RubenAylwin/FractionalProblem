#ifndef DISCRETE_MESH_TEST
#define DISCRETE_MESH_TEST
/*
 *TESTS FOR DISCRETE SPACES
 */
#include <Utilities.h>
#include <DiscreteSpaceMesh.h>
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
#include <Mesh.h>
#include <RegularP0Mesh.h>
#include <RegularP1Mesh.h>
namespace discreteMeshTests{

    static double tolerance = 1E-9;

    BOOST_AUTO_TEST_SUITE(DiscreteMesh)

    /*
     *TEST FOR MESH
     */

    BOOST_AUTO_TEST_SUITE(Mesh)
    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::MESH::ConstructorTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        msg(1) << "end Discrete::MESH::ConstructorTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ElementTest)
    {
        msg(1) << "start Discrete::MESH::ElementTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);

        for (int i = 0; i < 10; ++i){
            BOOST_CHECK_CLOSE(mesh.getElement(i).getA().getX(), i/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i).getA().getY(), 0.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i).getB().getX(), (i+1.)/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i).getB().getY(), 0.0, tolerance);
        }
        
        msg(1) << "end Discrete::MESH::ElementTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(ElementTest2)
    {
        msg(1) << "start Discrete::MESH::ElementTest2" << endMsg;
        StraightCurve curve(2.87);
        MeshCurve1D mesh(10, curve);

        for (int i = 0; i < 10; ++i){
            BOOST_CHECK_CLOSE(mesh.getElement(i).getA().getX(), 2.87*i/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i).getA().getY(), 0.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i).getB().getX(), 2.87*(i+1.)/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i).getB().getY(), 0.0, tolerance);
        }
        
        msg(1) << "end Discrete::MESH::ElementTest2" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(ElementEvaluationTest)
    {
        msg(1) << "start Discrete::MESH::ElementEvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);

        for (int i = 0; i < 10; ++i){
            BOOST_CHECK_CLOSE(mesh.getElement(i)(0.0).getX(), i/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i)(0.0).getY(), 0.0, tolerance);

            BOOST_CHECK_CLOSE(mesh.getElement(i)(0.5).getX(), (i+0.5)/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i)(0.5).getY(), 0.0, tolerance);
            
            BOOST_CHECK_CLOSE(mesh.getElement(i)(1.0).getX(), (i+1)/10.0, tolerance);
            BOOST_CHECK_CLOSE(mesh.getElement(i)(1.0).getY(), 0.0, tolerance);

        }
        
        msg(1) << "end Discrete::MESH::ElementEvaluationTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()
    
    /*
     *TEST FOR P0 SPACE
     */
    BOOST_AUTO_TEST_SUITE(P0)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P0::ConstructorTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP0Mesh_1D regular(mesh);
        msg(1) << "end Discrete::P0::ConstructorTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P0::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP0Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(9);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9,0.95).real(), tolerance);
        msg(1) << "end Discrete::P0::EvaluationTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationMeshTest)
    {
        msg(1) << "start Discrete::P0::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP0Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.15).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun.evaluate(2, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun.evaluate(3, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(9);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9, 0.95).real(), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun.evaluate(6, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun.evaluate(6, 0.95).real(), tolerance);

        msg(1) << "end Discrete::P0::EvaluationTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P0::InvalidFunctionTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP0Mesh_1D space(mesh);

        BOOST_REQUIRE_THROW(space.basisFunction(10), std::invalid_argument);
        msg(1) << "end Discrete::P0::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P0::TestingConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return 1. + 0*s;});
        auto testingResult = space.testAgainstBasis(function);
        
        for (const auto val : testingResult) {
            BOOST_CHECK_CLOSE(val.real(), 1.0/20, tolerance);
        }
        msg(1) << "end Discrete::P0::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingPiecewiseConstantFunction)
    {
        msg(1) << "start Discrete::P0::TestingPiecewiseConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return t > 0.5 ?  10 : 1;});
        auto testingResult = space.testAgainstBasis(function);
        
        for (unsigned i = 0; i < 10; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), 1.0/20, tolerance);
        }

        for (unsigned i = 10; i < 20; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), 0.5, tolerance);
        }

        msg(1) << "end Discrete::P0::TestingPiecewiseConstantFunction" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunction1)
    {
        msg(1) << "start Discrete::P0::DiscreteFunction1" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(20, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 1.0, tolerance);
        //BEM::plotFunction("testing", *fun, mesh);
        msg(1) << "end Discrete::P0::DiscreteFunction1" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunction2)
    {
        msg(1) << "start Discrete::P0::DiscreteFunction2" << endMsg;
        TrigonometricCurve curve(1.0, 0, {0.1}, {-0.2});
        MeshCurve1D mesh(20, curve);
        RegularP0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(20, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 1.0, tolerance);
        //BEM::plotFunction("testing", *fun, mesh);
        msg(1) << "end Discrete::P0::DiscreteFunction2" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()

    /*
     *TEST FOR P0 SPACE
     */
    BOOST_AUTO_TEST_SUITE(P0Graded)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P0Graded::ConstructorTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 1./0.4);
        RegularP0Mesh_1D regular(mesh);
        msg(1) << "end Discrete::P0Graded::ConstructorTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P0Graded::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 1./0.5);
        RegularP0Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(9);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9,0.95).real(), tolerance);
        msg(1) << "end Discrete::P0Graded::EvaluationTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationMeshTest)
    {
        msg(1) << "start Discrete::P0Graded::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 0.4);
        RegularP0Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, firstFun.evaluate(0, 0.15).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun.evaluate(2, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0, firstFun.evaluate(3, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(9);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1, lastFun.evaluate(9, 0.95).real(), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun.evaluate(6, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0, lastFun.evaluate(6, 0.95).real(), tolerance);

        msg(1) << "end Discrete::P0::EvaluationTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P0Graded::InvalidFunctionTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 0.5);
        RegularP0Mesh_1D space(mesh);

        BOOST_REQUIRE_THROW(space.basisFunction(10), std::invalid_argument);
        msg(1) << "end Discrete::P0Graded::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P0Graded::TestingConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(20, curve, 0.4);
        RegularP0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return 1. + 0*s;});
        auto testingResult = space.testAgainstBasis(function);
        
        for (int i = 0; i< testingResult.size(); ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), std::pow((i+1.)/20., 0.4) - std::pow((i)/20., 0.4), tolerance);
        }
        msg(1) << "end Discrete::P0Graded::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingPiecewiseConstantFunction)
    {
        msg(1) << "start Discrete::P0Graded::TestingPiecewiseConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(20, curve, 0.5);
        RegularP0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return t > std::pow(10./20., 0.5) ?  10 : 1;});
        auto testingResult = space.testAgainstBasis(function);
        
        for (unsigned i = 0; i < 10; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), (std::pow((i+1.)/20., 0.5) - std::pow((i)/20., 0.5)), tolerance);
        }

        for (unsigned i = 10; i < 20; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), 10*(std::pow((i+1.)/20., 0.5) - std::pow((i)/20., 0.5)), tolerance);
        }

        msg(1) << "end Discrete::P0Graded::TestingPiecewiseConstantFunction" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunction1)
    {
        msg(1) << "start Discrete::P0Graded::DiscreteFunction1" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(20, curve, 0.4);
        RegularP0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(20, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 1.0, tolerance);
        msg(1) << "end Discrete::P0Graded::DiscreteFunction1" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunction2)
    {
        msg(1) << "start Discrete::P0Graded::DiscreteFunction2" << endMsg;
        TrigonometricCurve curve(1.0, 0, {0.1}, {-0.2});
        MeshCurveGraded1D mesh(20, curve, 0.7);
        RegularP0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(20, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 1.0, tolerance);
        //BEM::plotFunction("testing", *fun, mesh);
        msg(1) << "end Discrete::P0Graded::DiscreteFunction2" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()


    /*
     *TEST FOR P1 SPACE
     */
    BOOST_AUTO_TEST_SUITE(P1_0)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P1_0::ConstructorTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP1_0Mesh_1D regular(mesh);
        msg(1) << "end Discrete::P1_0::ConstructorTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P1_0::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP1_0Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(0.05, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.15, firstFun.evaluate(0, 0.15).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, firstFun.evaluate(1, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.85, firstFun.evaluate(1, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(8);
        BOOST_CHECK_CLOSE(0.05, lastFun.evaluate(8, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, lastFun.evaluate(8,0.95).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.05, lastFun.evaluate(9,0.95).real(), tolerance);

        msg(1) << "end Discrete::P1_0::EvaluationTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P1_0::InvalidFunctionTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP1_0Mesh_1D space(mesh);

        BOOST_REQUIRE_THROW(space.basisFunction(9), std::invalid_argument);
        msg(1) << "end Discrete::P1_0::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P1_0::TestingConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return 1. + 0*s;});
        auto testingResult = space.testAgainstBasis(function);
        
        for (const auto val : testingResult) {
            BOOST_CHECK_CLOSE(val.real(), 1.0/20, tolerance);
        }
        msg(1) << "end Discrete::P1_0::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingPolyFunction)
    {
        msg(1) << "start Discrete::P1_0::TestingPolyFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return t*t;});
        auto testingResult = space.testAgainstBasis(function);
        for (int i = 0; i < 19; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), (12.*(i+1)*(i+1) + 2.)/96000., tolerance);
        }
        
        msg(1) << "end Discrete::P1_0::TestingPolyFunction" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunction1)
    {
        msg(1) << "start Discrete::P1_0::DiscreteFunction1" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 0.25, tolerance);
        msg(1) << "end Discrete::P1_0::DiscreteFunction1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunction2)
    {
        msg(1) << "start Discrete::P1_0::DiscreteFunction2" << endMsg;
        TrigonometricCurve curve(1.0, 0, {0.1}, {-0.2});
        MeshCurve1D mesh(20, curve);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), .25, tolerance);
        msg(1) << "end Discrete::P1_0::DiscreteFunction2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunctionDerivative)
    {
        msg(1) << "start Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[1]=1.0+BEM::I*2.0;
        base[10]=-1.0+BEM::I*2.0;
        base[11]=-1.0-BEM::I*2.0;
        base[12]=3.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        auto &dFun = fun->derivative();
        BEM::plotFunction("fun", *fun, mesh);
        BEM::plotFunction("der", dFun, mesh);
        msg(1) << "end Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunctionFracDer)
    {
        msg(1) << "start Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[0]=1.0+BEM::I*2.0;
        auto baseFun = space.generateFunction(base);
        std::vector<BEM::Complex> base2(19, 0.0);
        base2[17]=1.0+BEM::I*2.0;
        auto baseFun2 = space.generateFunction(base2);
        DiscreteFunctionMesh *fun = baseFun.get();
        DiscreteFunctionMesh *fun2 = baseFun2.get();
        RegularP1_0Mesh_1D::P1Function *pFun = dynamic_cast<RegularP1_0Mesh_1D::P1Function*>(fun);
        RegularP1_0Mesh_1D::P1Function *pFun2 = dynamic_cast<RegularP1_0Mesh_1D::P1Function*>(fun2);
        auto &dFunL = pFun->derivative(70);
        auto &dFunR = pFun2->derivative(-70);
        // BEM::plotFunction("fun", *fun, mesh);
        BEM::plotFunction("derL0", dFunL, mesh);
        BEM::plotFunction("derR0", dFunR, mesh);
        msg(1) << "end Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
    }

    BOOST_AUTO_TEST_SUITE_END()



    /*
     *TEST FOR P1 SPACE
     */
    BOOST_AUTO_TEST_SUITE(P1)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P1::ConstructorTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP1_Mesh_1D regular(mesh);
        msg(1) << "end Discrete::P1::ConstructorTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P1::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP1_Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(1.-0.05, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(1.-0.15, firstFun.evaluate(0, 0.15).real(), tolerance);
        BOOST_CHECK_CLOSE(0., firstFun.evaluate(1, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0., firstFun.evaluate(1, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(10);
        BOOST_CHECK_CLOSE(0., lastFun.evaluate(8, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0., lastFun.evaluate(8,0.95).real(), tolerance);
        BOOST_CHECK_CLOSE(0.05, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, lastFun.evaluate(9,0.95).real(), tolerance);

        msg(1) << "end Discrete::P1::EvaluationTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P1::InvalidFunctionTest" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(10, curve);
        RegularP1_Mesh_1D space(mesh);

        BOOST_REQUIRE_THROW(space.basisFunction(11), std::invalid_argument);
        msg(1) << "end Discrete::P1::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P1::TestingConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return 1. + 0*s;});
        auto testingResult = space.testAgainstBasis(function);


        for (int i = 0; i < 21; ++i) {
            auto val = testingResult[i]*((i == 0 or i == 20) ? 2. : 1.);
            BOOST_CHECK_CLOSE(val.real(), 1.0/20, tolerance);
        }
        msg(1) << "end Discrete::P1::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingPolyFunction)
    {
        msg(1) << "start Discrete::P1::TestingPolyFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return t*t;});
        auto testingResult = space.testAgainstBasis(function);
        BOOST_CHECK_CLOSE(testingResult[0].real(), std::pow(1./20., 3.)/3.-std::pow(1./20., 3.)/4., tolerance);
        BOOST_CHECK_CLOSE(testingResult[20].real(), (1.-std::pow(19./20., 4.))*20./4. - 19*(1.-std::pow(19./20., 3.))/3., tolerance);
        for (int i = 1; i < 20; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), (12.*(i)*(i) + 2.)/96000., tolerance);
        }
        
        msg(1) << "end Discrete::P1::TestingPolyFunction" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunction1)
    {
        msg(1) << "start Discrete::P1::DiscreteFunction1" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(21, 0.0);
        base[2]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 0.25, tolerance);
        msg(1) << "end Discrete::P1::DiscreteFunction1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunction2)
    {
        msg(1) << "start Discrete::P1::DiscreteFunction2" << endMsg;
        TrigonometricCurve curve(1.0, 0, {0.1}, {-0.2});
        MeshCurve1D mesh(20, curve);
        RegularP1_Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(21, 0.0);
        base[2]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), .25, tolerance);
        msg(1) << "end Discrete::P1::DiscreteFunction2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunctionDerivative)
    {
        msg(1) << "start Discrete::P1::DiscreteFunctionDerivative" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(21, 0.0);
        base[0]=1.0+BEM::I*2.0;
        base[10]=-1.0+BEM::I*2.0;
        base[11]=-1.0-BEM::I*2.0;
        base[12]=3.0+BEM::I*2.0;
        base[20]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        auto &dFun = fun->derivative();
        BEM::plotFunction("fun", *fun, mesh);
        BEM::plotFunction("der", dFun, mesh);
        msg(1) << "end Discrete::P1::DiscreteFunctionDerivative" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunctionFracDer)
    {
        msg(1) << "start Discrete::P1::DiscreteFunctionDerivative" << endMsg;
        StraightCurve curve(1);
        MeshCurve1D mesh(20, curve);
        RegularP1_Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(21, 0.0);
        base[0]=1.0+BEM::I*2.0;
        auto baseFun = space.generateFunction(base);
        std::vector<BEM::Complex> base2(21, 0.0);
        base2[20]=1.0+BEM::I*2.0;
        auto baseFun2 = space.generateFunction(base2);
        DiscreteFunctionMesh *fun = baseFun.get();
        DiscreteFunctionMesh *fun2 = baseFun2.get();
        RegularP1_Mesh_1D::P1Function *pFun = dynamic_cast<RegularP1_Mesh_1D::P1Function*>(fun);
        RegularP1_Mesh_1D::P1Function *pFun2 = dynamic_cast<RegularP1_Mesh_1D::P1Function*>(fun2);
        auto &dFunL = pFun->derivative(70);
        auto &dFunR = pFun2->derivative(-70);
        BEM::plotFunction("fun", *fun, mesh);
        BEM::plotFunction("derL", dFunL, mesh);
        BEM::plotFunction("derR", dFunR, mesh);
        msg(1) << "end Discrete::P1::DiscreteFunctionDerivative" << endMsg;
    }
    
    BOOST_AUTO_TEST_SUITE_END()
    
    /*
     *TEST FOR P1 SPACE
     */
    BOOST_AUTO_TEST_SUITE(P1Graded)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        msg(1) << "start Discrete::P1Graded::ConstructorTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 0.8);
        RegularP1_0Mesh_1D regular(mesh);
        msg(1) << "end Discrete::P1Graded::ConstructorTest" << endMsg;
    }

    
    BOOST_AUTO_TEST_CASE(EvaluationTest)
    {
        msg(1) << "start Discrete::P1Graded::EvaluationTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 0.7);
        RegularP1_0Mesh_1D space(mesh);

        auto &firstFun = space.basisFunction(0);
        BOOST_CHECK_CLOSE(0.05, firstFun.evaluate(0, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.15, firstFun.evaluate(0, 0.15).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, firstFun.evaluate(1, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.85, firstFun.evaluate(1, 0.15).real(), tolerance);

        auto &lastFun = space.basisFunction(8);
        BOOST_CHECK_CLOSE(0.05, lastFun.evaluate(8, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, lastFun.evaluate(8,0.95).real(), tolerance);
        BOOST_CHECK_CLOSE(0.95, lastFun.evaluate(9, 0.05).real(), tolerance);
        BOOST_CHECK_CLOSE(0.05, lastFun.evaluate(9,0.95).real(), tolerance);

        msg(1) << "end Discrete::P1Graded::EvaluationTest" << endMsg;
    }
    
    BOOST_AUTO_TEST_CASE(InvalidFunctionTest)
    {
        msg(1) << "start Discrete::P1_0::InvalidFunctionTest" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(10, curve, 0.9);
        RegularP1_0Mesh_1D space(mesh);

        BOOST_REQUIRE_THROW(space.basisFunction(9), std::invalid_argument);
        msg(1) << "end Discrete::P1_0::InvalidFunctionTest" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(TestingConstantFunction)
    {
        msg(1) << "start Discrete::P1_0::TestingConstantFunction" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(20, curve, 0.8);
        RegularP1_0Mesh_1D space(mesh);

        ExplicitScalarFunction_2D function([&](double t, double s){ return 1. + 0*s;});
        auto testingResult = space.testAgainstBasis(function);
        
        for (int i = 0; i < 19; ++i) {
            BOOST_CHECK_CLOSE(testingResult[i].real(), 0.5*(std::pow((i+2.)/20., 0.8) - std::pow(i/20., 0.8)), tolerance);
        }
        msg(1) << "end Discrete::P1_0::TestingConstantFunction" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunction1)
    {
        msg(1) << "start Discrete::P1Graded::DiscreteFunction1" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(20, curve, 0.8);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), 0.25, tolerance);
        msg(1) << "end Discrete::P1Graded::DiscreteFunction1" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunction2)
    {
        msg(1) << "start Discrete::P1Graded::DiscreteFunction2" << endMsg;
        TrigonometricCurve curve(1.0, 0, {0.1}, {-0.2});
        MeshCurveGraded1D mesh(20, curve, 0.7);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[1]=1.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        BOOST_CHECK_CLOSE(fun->evaluate(0, 0.5).real(), 0.0, tolerance);
        BOOST_CHECK_CLOSE(fun->evaluate(1, 0.25).real(), .25, tolerance);
        msg(1) << "end Discrete::P1Graded::DiscreteFunction2" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(DiscreteFunctionDerivative)
    {
        msg(1) << "start Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
        StraightCurve curve(1);
        MeshCurveGraded1D mesh(20, curve, .2);
        RegularP1_0Mesh_1D space(mesh);
        std::vector<BEM::Complex> base(19, 0.0);
        base[1]=1.0+BEM::I*2.0;
        base[10]=-1.0+BEM::I*2.0;
        base[11]=-1.0-BEM::I*2.0;
        base[12]=3.0+BEM::I*2.0;
        auto fun = space.generateFunction(base);
        auto &dFun = fun->derivative();
        BEM::plotFunction("fun", *fun, mesh);
        BEM::plotFunction("der", dFun, mesh);
        msg(1) << "end Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(DiscreteFunctionDerivativeProj)
    {
        msg(1) << "start Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
        StraightCurve curve(1);
        // MeshCurve1D mesh(200, curve);
        MeshCurveGraded1D mesh(200, curve, 1.6);
        RegularP1_0Mesh_1D space(mesh);
        auto fun = space.generateFunction(ExplicitScalarFunction_2D([](double t, double s){ return t*(1.-t);}));
        auto &dFun = fun->derivative();
        BEM::plotFunction("fun1", *fun, mesh);
        BEM::plotFunction("der1", dFun, mesh);
        msg(1) << "end Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
    }


    // BOOST_AUTO_TEST_CASE(DiscreteFunctionFracDer)
    // {
    //     msg(1) << "start Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
    //     StraightCurve curve(1);
    //     MeshCurve1D mesh(20, curve);
    //     RegularP1_0Mesh_1D space(mesh);
    //     std::vector<BEM::Complex> base(19, 0.0);
    //     base[1]=1.0+BEM::I*2.0;
    //     auto baseFun = space.generateFunction(base);
    //     std::vector<BEM::Complex> base2(19, 0.0);
    //     base2[17]=1.0+BEM::I*2.0;
    //     auto baseFun2 = space.generateFunction(base2);
    //     DiscreteFunctionMesh *fun = baseFun.get();
    //     DiscreteFunctionMesh *fun2 = baseFun2.get();
    //     RegularP1_0Mesh_1D::P1Function *pFun = dynamic_cast<RegularP1_0Mesh_1D::P1Function*>(fun);
    //     RegularP1_0Mesh_1D::P1Function *pFun2 = dynamic_cast<RegularP1_0Mesh_1D::P1Function*>(fun2);
    //     auto &dFunL = pFun->derivative(70);
    //     auto &dFunR = pFun2->derivative(-70);
    //     // BEM::plotFunction("fun", *fun, mesh);
    //     BEM::plotFunction("derL", dFunL, mesh);
    //     BEM::plotFunction("derR", dFunR, mesh);
    //     msg(1) << "end Discrete::P1_0::DiscreteFunctionDerivative" << endMsg;
    // }

    BOOST_AUTO_TEST_SUITE_END()


    
    BOOST_AUTO_TEST_SUITE_END()
}
#endif
