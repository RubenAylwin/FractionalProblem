#ifndef RL_TEST
#define RL_TEST
#include <FractionalParabolic.h>
#include <cmath>
#include <RiemannLiouville.h>
#include <SobolevSlobodeckii.h>
#include <ScalarValuedFunction.h>
#include <ParametrizedFunction.h>
#include <RiemannLiouvilleProblem.h>
#include <L2Operator.h>
#include <DiscreteSpace.h>
#include <Curve.h>
#include <vector>
#include <MyTypes.h>
#include <Utilities.h>
#include <RegularP1Mesh.h>
#include <RegularP0Mesh.h>
#include <algorithm>
#include <random>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <random>

namespace riemannLiouvilleTests {
    static double tolerance = 1e-3;
    BOOST_AUTO_TEST_SUITE(Riemann)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
        ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
        unsigned size = std::pow(2, 3);
        int order = 85;
        MeshCurve1D mesh(size, curve);
        RegularP1_0Mesh_1D space(mesh);
        RiemannLiouvilleMesh rlm(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    }

    BOOST_AUTO_TEST_CASE(ConstructorTest2)
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        PolynomialFunction_1D d(BEM::CVector{1.0});
        TrigonometricFunction_1D q(1, BEM::CVector{0.0, 0.0});
        unsigned size = std::pow(2, 3);
        int order = 85;
        MeshCurve1D mesh(size, curve);
        RegularP1_0Mesh_1D space(mesh);
        RiemannLiouvilleMesh rlm(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    }

    BOOST_AUTO_TEST_CASE(DifferentFunTypeTest)
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
        ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return t;});
        
        unsigned size = std::pow(2, 10);
        int order = 85;
        MeshCurve1D mesh(size, curve);
        RegularP1_0Mesh_1D space(mesh);
        RiemannLiouvilleMesh rlm1(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);

        TrigonometricFunction_1D d2(1, BEM::CVector{1., 0.0});
        PolynomialFunction_1D q2(BEM::CVector{0., 1.0});

        RiemannLiouvilleMesh rlm2(space, RiemannLiouvilleMesh::Side::LEFT, order, d2, q2);

        rlm1.assembleMassMatrix();
        rlm2.assembleMassMatrix();

        BOOST_CHECK_CLOSE((rlm1.getMatrix() - rlm2.getMatrix()).operatorNorm(), 0.0, tolerance);
    }

    // BOOST_AUTO_TEST_CASE(ProblemMeshCheck, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

    //     int order = 85;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;
    //     ExplicitScalarFunction_2D analyticSolution([=](double t, double s) {return (-std::pow(t, rOrder) + std::pow(t, rOrder-1))/(rOrder*std::tgamma(rOrder));});
    //     ExplicitScalarFunction_2D analyticSolutionDer([=](double t, double s) {return (-std::tgamma(rOrder + 1)/std::tgamma(rOrder+1-fOrder)*std::pow(t, fOrder) + std::tgamma(rOrder)/std::tgamma(fOrder)*std::pow(t, fOrder-1))/(rOrder*std::tgamma(rOrder));});

    //     unsigned size = std::pow(2, 6);
    //     MeshCurve1D mesh(size, curve);
    //     RegularP1_0Mesh_1D space(mesh);
    //     RiemannLiouvilleMesh rlm(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
        
    //     rlm.assembleMassMatrix();
    //     auto vec = space.testAgainstBasis(rhs);
    //     BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //     std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //     for (int i = 0; i < size-1; ++i) {
    //         solutionVec[i] = solutionVec1[i];
    //     }
    //     auto solution = space.generateFunction(solutionVec);
    //     auto &solutionDerL = solution->derivative(order);
    //     //auto &solutionDerR = solution->derivative(-order);
    //     BEM::plotFunction("solutionCheck", *solution, mesh);
    //     BEM::plotFunction("solutionDerLeftCheck", solutionDerL, mesh);

    //     RiemannLiouvilleProblem RLP(order, space, q, d, rhs);
    //     RLP.buildDiscrete();
    //     RLP.solve();
    //     std::vector<BEM::Complex> solutionVecP(size-1, 0.0);
    //     for (int i = 0; i < size-1; ++i) {
    //         solutionVecP[i] = RLP.getSolutionVec()[i];
    //     }
    //     auto solutionP = space.generateFunction(solutionVecP);
    //     auto &solutionDerLP = solutionP->derivative(order);
    //     BEM::plotFunction("solutionCheckProb", *solutionP, mesh);
    //     BEM::plotFunction("solutionDerLeftCheckProb", solutionDerLP, mesh);
        
    // }


    // // BOOST_AUTO_TEST_CASE(ProblemMeshCheckCaputo, *boost::unit_test::disabled())
    // // {
    // //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    // //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    // //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    // //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 0. - 0.0;});

    // //     int order = 95;
    // //     double fOrder = order/100.;
    // //     double rOrder = 2*fOrder;
    // //     ExplicitScalarFunction_2D analyticSolution([=](double t, double s) {return (-std::pow(t, rOrder) + std::pow(t, rOrder-1))/(rOrder*std::tgamma(rOrder));});
    // //     ExplicitScalarFunction_2D analyticSolutionDer([=](double t, double s) {return (-std::tgamma(rOrder + 1)/std::tgamma(rOrder+1-fOrder)*std::pow(t, fOrder) + std::tgamma(rOrder)/std::tgamma(fOrder)*std::pow(t, fOrder-1))/(rOrder*std::tgamma(rOrder));});

    // //     unsigned size = std::pow(2, 6);
    // //     MeshCurve1D mesh(size, curve);
    // //     RegularP1_Mesh_1D space(mesh);
    // //     RiemannLiouvilleMesh rlm(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
        
    // //     rlm.assembleMassMatrix();
    // //     auto vec = space.testAgainstBasis(rhs);
    // //     vec[1]=-1.0;
    // //     vec[size]=1.0;
    // //     vec[size-1]=-1.0;
    // //     BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    // //     std::vector<BEM::Complex> solutionVec(size+1, 0.0);
    // //     for (int i = 0; i < size+1; ++i) {
    // //         solutionVec[i] = solutionVec1[i];
    // //     }
    // //     auto solution = space.generateFunction(solutionVec);
    // //     auto &solutionDerL = solution->derivative(order);
    // //     BEM::plotFunction("solutionCheck", *solution, mesh);
    // //     BEM::plotFunction("solutionDerLeftCheck", solutionDerL, mesh);
    // // }


    // BOOST_AUTO_TEST_CASE(FactoryMeshCheck, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

    //     int order = 85;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;

    //     unsigned size = std::pow(2, 7);
    //     MeshCurve1D mesh(size, curve);
    //     RegularP1_0Mesh_1D space(mesh);

    //     RiemannLiouvilleProblem RLP(order, space, q, d, rhs);
    //     RLP.buildDiscrete();
    //     RLP.solve();
    //     std::vector<BEM::Complex> solutionVecP(size-1, 0.0);
    //     for (int i = 0; i < size-1; ++i) {
    //         solutionVecP[i] = RLP.getSolutionVec()[i];
    //     }
    //     auto solutionP = space.generateFunction(solutionVecP);
    //     auto &solutionDerLP = solutionP->derivative(order);
    //     BEM::plotFunction("solutionCheckProb", *solutionP, mesh);
    //     BEM::plotFunction("solutionDerLeftCheckProb", solutionDerLP, mesh);

    //     ExplicitScalarFunction_1D d1([](double t)->BEM::Complex {return t < 0.5 ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D q1([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_2D rhs1([](double t, double s)->BEM::Complex {return .5;});
    //     ExplicitScalarFunction_1D d2([](double t)->BEM::Complex {return t < 0.5 ? 0.0 : 1.0;});
    //     ExplicitScalarFunction_1D q2([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_2D rhs2([](double t, double s)->BEM::Complex {return .25;});
    //     RiemannLiouvilleMeshFactory factory(size, order, VectorFun_1D{q1, q2}, VectorFun_1D{d1, d2}, VectorFun_2D{rhs1, rhs2});
    //     auto factProb = factory.generateNewProblem(std::vector<double>{1., -1., 1., 1., 1., 2.});
    //     factProb->buildDiscrete();
    //     factProb->solve();
    //     BOOST_CHECK_CLOSE((RLP.getSolutionVec() - factProb->getSolutionVec()).norm(), 0.0, tolerance);
    // }

    // BOOST_AUTO_TEST_CASE(FactoryGreedyCheck, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

    //     int order = 90;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;

    //     unsigned size = std::pow(2, 8);
    //     ExplicitScalarFunction_1D d1([](double t)->BEM::Complex {return t < 0.25 ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D d2([](double t)->BEM::Complex {return (t >= 0.25 and t < 0.5) ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D d3([](double t)->BEM::Complex {return (t >= 0.5 and t < 0.75) ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D d4([](double t)->BEM::Complex {return t >= 0.75 ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D q1([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_2D rhs1([](double t, double s)->BEM::Complex {return t*(1.-t);});
    //     ExplicitScalarFunction_1D q2([](double t)->BEM::Complex {return t;});
    //     ExplicitScalarFunction_2D rhs2([](double t, double s)->BEM::Complex {return 1.;});
    //     RiemannLiouvilleMeshFactory factory(size, order, VectorFun_1D{}, VectorFun_1D{d1, d2, d3, d4}, VectorFun_2D{rhs1, rhs2});


    //     std::vector<BEM::Interval1D> limits{BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.1, 1), BEM::Interval1D(0.1, 1)};
    //     std::vector<std::vector<double>> pointsForTesting{};
    //     std::uniform_real_distribution<double> _unif(0., 1.);
    //     std::mt19937 _rng(1993);
    //     for(int i = 0; i < 200; ++i) {
    //         pointsForTesting.push_back(std::vector<double>{1+0.5*_unif(_rng), 1+0.5*_unif(_rng), 1+0.5*_unif(_rng), 1+0.5*_unif(_rng), .1+0.9*_unif(_rng), .1+0.9*_unif(_rng)});
    //     }
    //     auto est = [fOrder](std::vector<double> point) -> double {
    //         double min = *std::min_element(point.begin(), point.begin() + 4);
    //         double max = *std::max_element(point.begin(), point.begin() + 4);
            
    //         return 0.5*((max + min)*std::abs(std::cos(M_PI*fOrder)) - (max - min));};
    //     startTimer(training);
    //     factory.trainGreedy(limits, 10, 1e-4, est);
    //     stopTimer(training, "Training");
    //     startTimer(greedySolve);
    //     auto greedySolution = factory.greedySolve(std::vector<double>{1.0, 0.7, 1.3, 1.4, 1.0, 0.0});
    //     stopTimer(greedySolve, "Greedy solution");
    //     auto prob = factory.generateNewProblem(std::vector<double>{1.0, 0.7, 1.3, 1.4, 1.0, 0.0});
    //     startTimer(assemble);
    //     prob->buildDiscrete();
    //     stopTimer(assemble, "Assemble HFM");
    //     startTimer(solveFull);
    //     prob->solve();
    //     stopTimer(solveFull, "solution HFM");
    //     auto solution = prob->getSolutionVec();
        
    // }

    BOOST_AUTO_TEST_CASE(FactoryGreedyCheck33, *boost::unit_test::disabled())
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
        ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
        ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

        int order = 90;
        double fOrder = order/100.;
        double rOrder = 2*fOrder;

        unsigned size = std::pow(2, 8);
        ExplicitScalarFunction_1D d1([](double t)->BEM::Complex {return t < 0.25 ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D d2([](double t)->BEM::Complex {return (t >= 0.25 and t < 0.5) ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D d3([](double t)->BEM::Complex {return (t >= 0.5 and t < 0.75) ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D d4([](double t)->BEM::Complex {return t >= 0.75 ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D q1([](double t)->BEM::Complex {return 1.0;});
        ExplicitScalarFunction_2D rhs1([](double t, double s)->BEM::Complex {return t*(1.-t);});
        ExplicitScalarFunction_1D q2([](double t)->BEM::Complex {return t;});
        ExplicitScalarFunction_2D rhs2([](double t, double s)->BEM::Complex {return 1.;});
        RiemannLiouvilleMeshFactory factory(size, order, VectorFun_1D{}, VectorFun_1D{d1, d2, d3, d4}, VectorFun_2D{rhs1});


        std::vector<BEM::Interval1D> limits{BEM::Interval1D(1., 1.), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(1., 1)};
        auto est = [fOrder](std::vector<double> point) -> double {
            double min = *std::min_element(point.begin(), point.begin() + 4);
            double max = *std::max_element(point.begin(), point.begin() + 4);
            
            return 0.5*((max + min)*std::abs(std::cos(M_PI*fOrder)) - (max - min));};
        startTimer(training);
        factory.trainGreedy(limits, 10, 1e-4, est);
        stopTimer(training, "Training");
        // startTimer(greedySolve);
        // auto greedySolution = factory.greedySolve(std::vector<double>{1.0, 0.7, 1.3, 1.4, 1.0, 0.0});
        // stopTimer(greedySolve, "Greedy solution");
        // auto prob = factory.generateNewProblem(std::vector<double>{1.0, 0.7, 1.3, 1.4, 1.0, 0.0});
        // startTimer(assemble);
        // prob->buildDiscrete();
        // stopTimer(assemble, "Assemble HFM");
        // startTimer(solveFull);
        // prob->solve();
        // stopTimer(solveFull, "solution HFM");
        // auto solution = prob->getSolutionVec();
        
    }
    
    BOOST_AUTO_TEST_CASE(FactoryGreedyCheck2, *boost::unit_test::disabled())
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
        ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
        ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

        int order = BEM::getEnv<int>("TEST2_FORDER");
        order = order > 50 ? order : 60;
        double fOrder = order/100.;
        double rOrder = 2*fOrder;

        unsigned size = std::pow(2, 6);
        ExplicitScalarFunction_1D d1([](double t)->BEM::Complex {return t < 0.25 ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D d2([](double t)->BEM::Complex {return (t >= 0.25 and t < 0.5) ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D d3([](double t)->BEM::Complex {return (t >= 0.5 and t < 0.75) ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D d4([](double t)->BEM::Complex {return t >= 0.75 ? 1.0 : 0.0;});
        ExplicitScalarFunction_1D q1([](double t)->BEM::Complex {return 1.0;});
        ExplicitScalarFunction_2D rhs1([](double t, double s)->BEM::Complex {return t*(1.-t);});
        ExplicitScalarFunction_1D q2([](double t)->BEM::Complex {return t;});
        ExplicitScalarFunction_1D q3([](double t)->BEM::Complex {return (1.-t);});
        ExplicitScalarFunction_1D q4([](double t)->BEM::Complex {return t*t;});
        ExplicitScalarFunction_1D q5([](double t)->BEM::Complex {return (1.-t)*(1.-t);});
        ExplicitScalarFunction_2D rhs2([](double t, double s)->BEM::Complex {return 1.;});
        RiemannLiouvilleMeshFactory factory(size, order, VectorFun_1D{d}, VectorFun_1D{d}, VectorFun_2D{rhs1});
        std::vector<std::vector<double>> pointsForTesting{};
        std::uniform_real_distribution<double> _unif(0., 1.);
        std::mt19937 _rng(1993);
        auto est = [fOrder](std::vector<double> point) -> double {
            // double min = *std::min_element(point.begin() + 2, point.begin() + 4 + 2);
            // double max = *std::max_element(point.begin() + 2, point.begin() + 4 + 2);
            return std::abs(std::cos(M_PI*fOrder));
            // return 0.5*((max + min)*std::abs(std::cos(M_PI*fOrder)) - (max - min));
        };

        // std::vector<BEM::Interval1D> limits{BEM::Interval1D(0., 1.), BEM::Interval1D(0., 1.), BEM::Interval1D(1, 1), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(0.7, 1.3), BEM::Interval1D(1., 1.)};
        std::vector<BEM::Interval1D> limits{BEM::Interval1D(0., 5.), BEM::Interval1D(1., 1.),  BEM::Interval1D(1., 1.)};
        factory.trainGreedy(limits, 100, 1e-6, est);

    }
    

    // BOOST_AUTO_TEST_CASE(FactoryGreedyCheck3, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

    //     int order = 90;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;

    //     unsigned size = std::pow(2, 8);
    //     ExplicitScalarFunction_1D d1([](double t)->BEM::Complex {return t < 0.25 ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D d2([](double t)->BEM::Complex {return (t >= 0.25 and t < 0.5) ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D d3([](double t)->BEM::Complex {return (t >= 0.5 and t < 0.75) ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D d4([](double t)->BEM::Complex {return t >= 0.75 ? 1.0 : 0.0;});
    //     ExplicitScalarFunction_1D q1([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_2D rhs1([](double t, double s)->BEM::Complex {return t*(1.-t);});
    //     ExplicitScalarFunction_1D q2([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_2D rhs2([](double t, double s)->BEM::Complex {return 1.;});
    //     RiemannLiouvilleMeshFactory factory(size, order, VectorFun_1D{q1}, VectorFun_1D{q2}, VectorFun_2D{rhs2});


    //     std::vector<std::vector<double>> pointsForTesting{};
    //     std::uniform_real_distribution<double> _unif(0., 1.);
    //     std::mt19937 _rng(1993);
    //     for(int i = 0; i < 200; ++i) {
    //         pointsForTesting.push_back(std::vector<double>{1.+4.*_unif(_rng), 1., 1.});
    //     }
        
    //     std::vector<BEM::Interval1D> limits{BEM::Interval1D(1., 5), BEM::Interval1D(1, 1), BEM::Interval1D(1, 1)};
    //     auto est = [fOrder](std::vector<double> point) -> double {
    //         return std::abs(std::cos(M_PI*fOrder));
    //         // return 1;
    //     };
    //     factory.trainGreedy(limits, 100, 1e-6, est);

    // }
    
    // BOOST_AUTO_TEST_CASE(MatrixMesh, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return 1.0;});

    //     int order = 60;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;
    //     ExplicitScalarFunction_2D analyticSolution([=](double t, double s) {return (-std::pow(t, rOrder) + std::pow(t, rOrder-1))/(rOrder*std::tgamma(rOrder));});
    //     ExplicitScalarFunction_2D analyticSolutionDer([=](double t, double s) {return (-std::tgamma(rOrder + 1)/std::tgamma(rOrder+1-fOrder)*std::pow(t, fOrder) + std::tgamma(rOrder)/std::tgamma(fOrder)*std::pow(t, fOrder-1))/(rOrder*std::tgamma(rOrder));});
    //     for (int s = 1; s <= 10; ++s) {   
    //         unsigned size = std::pow(2, s);
    //         MeshCurve1D mesh(size, curve);
    //         RegularP1_0Mesh_1D space(mesh);
    //         RiemannLiouvilleMesh rlm(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //         rlm.assembleMassMatrix();
    //         auto &matrix = rlm.getMatrix();
    //         Eigen::BDCSVD<BEM::Matrix, Eigen::ComputeThinU> SVD(matrix);
    //         auto vec = space.testAgainstBasis(rhs);
    //         BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //         std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //         for (int i = 0; i < size-1; ++i) {
    //             solutionVec[i] = solutionVec1[i];
    //         }
    //         auto solution = space.generateFunction(solutionVec);
    //         auto &solutionDerL = solution->derivative(order);
    //         BEM::plotFunction("solution", *solution, mesh);
    //         BEM::plotFunction("solutionDerLeft", solutionDerL, mesh);
    //     }
    // }

    // BOOST_AUTO_TEST_CASE(MatrixMesh2, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return 0.0;});
    //     ExplicitScalarFunction_2D rhs([](double t, double s)->BEM::Complex {return t*(1.-t);});

    //     int order = 60;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;
    //     ExplicitScalarFunction_2D analyticSolution([=](double t, double s) {return (std::pow(t, rOrder-1) - std::pow(t, rOrder+1))/(std::tgamma(rOrder+2)) - 2.*(std::pow(t, rOrder-1) - std::pow(t, rOrder+2))/(std::tgamma(rOrder+3));});
    //     ExplicitScalarFunction_2D analyticSolutionDer([=](double t, double s) {return (std::pow(t, rOrder-fOrder-1)*std::tgamma(rOrder)/std::tgamma(rOrder-fOrder) - std::pow(t, rOrder-fOrder+1)*std::tgamma(rOrder+2)/std::tgamma(rOrder+2-fOrder))/(std::tgamma(rOrder+2)) -
    //                 2.*(std::pow(t, rOrder-fOrder-1)*std::tgamma(rOrder)/std::tgamma(rOrder-fOrder) - std::pow(t, rOrder-fOrder+2)*std::tgamma(rOrder+3)/std::tgamma(rOrder-fOrder+3))/(std::tgamma(rOrder+3));});

    //     for (int s = 1; s <= 10; ++s) {   
    //         unsigned size = std::pow(2, s);
    //         MeshCurve1D mesh(size, curve);
    //         RegularP1_0Mesh_1D space(mesh);
    //         RiemannLiouvilleMesh rlm(space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //         startTimer(assemble);
    //         rlm.assembleMassMatrix();
    //         auto vec = space.testAgainstBasis(rhs);
    //         stopTimer(assemble,"Assembled "+std::to_string(s));
    //         startTimer(solve)
    //         BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //         stopTimer(solve,"Solved "+std::to_string(s));
    //         std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //         for (int i = 0; i < size-1; ++i) {
    //             solutionVec[i] = solutionVec1[i];
    //         }
    //         auto solution = space.generateFunction(solutionVec);
    //         auto &solutionDer = solution->derivative(order);
    //     }
    // }


    // BOOST_AUTO_TEST_CASE(MatrixMesh3, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return t < 0.5 ? 5.0 : 3.0;});
    //     int order = 60;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;
    //     // ExplicitScalarFunction_1D analyticSolution([=](double t) {return (-std::pow(t, rOrder) + std::pow(t, rOrder-1))/(rOrder*std::tgamma(rOrder));});
    //     // // ExplicitScalarFunction_2D analyticSolutionDer([=](double t, double s) {return (-std::tgamma(rOrder + 1)/std::tgamma(rOrder+1-fOrder)*std::pow(t, fOrder) + std::tgamma(rOrder)/std::tgamma(fOrder)*std::pow(t, fOrder-1))/(rOrder*std::tgamma(rOrder));});
    //     // BEM::plotFunction("aSol", analyticSolution);
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {
    //         return t < 0.5 ? -2.0 : 8.0;});
    //     ExplicitScalarFunction_2D rhs([&](double t, double s)->BEM::Complex {return 1.;});
    //     int size = std::pow(2, 1);
    //     int power = 11;
    //     int sizeOK = std::pow(2, power);
    //     MeshCurveGraded1D meshOK = MeshCurveGraded1D(sizeOK, curve, 1./(rOrder-1.));
    //     // auto meshOK = MeshCurve1D(sizeOK, curve);
    //     // RegularP1_0Mesh_1D spaceOK(meshOK);
    //     // RegularP0Mesh_1D spaceOK0(meshOK);
    //     // RiemannLiouvilleMesh rlmOK(spaceOK, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //     // auto vecOK = spaceOK.testAgainstBasis(rhs);
    //     // BEM::ColVector solutionVecOK1(rlmOK.getMatrix().colPivHouseholderQr().solve(vecOK));
    //     // std::vector<BEM::Complex> solutionVecOK(sizeOK-1, 0.0);
    //     // for (int i = 0; i < sizeOK-1; ++i) {
    //     //     solutionVecOK[i] = solutionVecOK1[i];
    //     // }
    //     // auto solutionOK = spaceOK.generateFunction(solutionVecOK);
    //     // auto &solutionDerOKL = solutionOK->derivative(order);
    //     // auto &solutionDerOKR = solutionOK->derivative(-order);
    //     // MeshCurve1D *meshPrev = new MeshCurve1D(size, curve);
    //     // RegularP1_0Mesh_1D *spacePrev = new RegularP1_0Mesh_1D(*meshPrev);
    //     // RiemannLiouvilleMesh rlm(*spacePrev, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //     // auto vec = spacePrev->testAgainstBasis(rhs);
    //     // BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //     // std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //     // for (int i = 0; i < size-1; ++i) {
    //     //     solutionVec[i] = solutionVec1[i];
    //     // }
    //     // auto solutionPrev = spacePrev->generateFunction(solutionVec);
    //     // BEM::plotFunction("solutionOKP"+std::to_string(order-50), *solutionOK, meshOK);
    //     // BEM::plotFunction("solutionDerOKP", solutionDerOKL, meshOK);
    //     for (int s = 1; s <= power-1; ++s) {
    //         size = std::pow(2, s);
    //         MeshCurve1D *mesh = new MeshCurve1D(size, curve);
    //         RegularP1_0Mesh_1D *space = new RegularP1_0Mesh_1D(*mesh);
    //         RiemannLiouvilleMesh rlm(*space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //         auto &matrix = rlm.getMatrix();
    //         Eigen::BDCSVD<BEM::Matrix, Eigen::ComputeThinU> SVD(matrix);
    //         // auto vec = space->testAgainstBasis(rhs);
    //         // BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //         // std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //         // for (int i = 0; i < size-1; ++i) {
    //         //     solutionVec[i] = solutionVec1[i];
    //         // }
    //         // auto solution = space->generateFunction(solutionVec);
    //         // auto &solutionDerL = solution->derivative(order);
    //         // BEM::plotFunction("sol"+std::to_string(s), *solution, *mesh);
    //         // BEM::plotFunction("der"+std::to_string(s), solutionDerL, *mesh);
    //         // std::cout << solution->L2Error(ExplicitScalarFunction_2D([&](double t, double s){return (*solutionOK)(t);})) << " " <<
    //         //           solutionDerL.L2Error(ExplicitScalarFunction_2D([&](double t, double s){return solutionDerOKL(t);})) << std::endl;
    //         delete mesh;
    //         delete space;
    //         // delete meshPrev;
    //         // delete spacePrev;
    //         // meshPrev = mesh;
    //         // spacePrev = space;
    //     }
    // }

    // BOOST_AUTO_TEST_CASE(MatrixMesh4, *boost::unit_test::disabled())
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 4.+std::sin(2*M_PI*t);});
    //     int order = 90;
    //     double fOrder = order/100.;
    //     double rOrder = 2*fOrder;
        // ExplicitScalarFunction_1D analyticSolution([=](double t) {return (-std::pow(t, rOrder) + std::pow(t, rOrder-1))/(rOrder*std::tgamma(rOrder));});
        // // ExplicitScalarFunction_2D analyticSolutionDer([=](double t, double s) {return (-std::tgamma(rOrder + 1)/std::tgamma(rOrder+1-fOrder)*std::pow(t, fOrder) + std::tgamma(rOrder)/std::tgamma(fOrder)*std::pow(t, fOrder-1))/(rOrder*std::tgamma(rOrder));});
        // BEM::plotFunction("aSol", analyticSolution);
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {
    //         return std::cos(2*M_PI*t);});
    //     ExplicitScalarFunction_2D rhs([&](double t, double s)->BEM::Complex {return 1.;});
    //     int size = std::pow(2, 1);
    //     int power = 9;
    //     int sizeOK = std::pow(2, power);
    //     MeshCurveGraded1D meshOK = MeshCurveGraded1D(sizeOK, curve, 1./(rOrder-1.));
    //     // // auto meshOK = MeshCurve1D(sizeOK, curve);
    //     RegularP1_0Mesh_1D spaceOK(meshOK);
    //     // RegularP0Mesh_1D spaceOK0(meshOK);
    //     RiemannLiouvilleMesh rlmOK(spaceOK, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //     auto vecOK = spaceOK.testAgainstBasis(rhs);
    //     BEM::ColVector solutionVecOK1(rlmOK.getMatrix().colPivHouseholderQr().solve(vecOK));
    //     std::vector<BEM::Complex> solutionVecOK(sizeOK-1, 0.0);
    //     for (int i = 0; i < sizeOK-1; ++i) {
    //         solutionVecOK[i] = solutionVecOK1[i];
    //     }
    //     auto solutionOK = spaceOK.generateFunction(solutionVecOK);
    //     auto &solutionDerOKL = solutionOK->derivative(order);
    //     auto &solutionDerOKR = solutionOK->derivative(-order);
    //     MeshCurve1D *meshPrev = new MeshCurve1D(size, curve);
    //     RegularP1_0Mesh_1D *spacePrev = new RegularP1_0Mesh_1D(*meshPrev);
    //     RiemannLiouvilleMesh rlm(*spacePrev, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //     auto vec = spacePrev->testAgainstBasis(rhs);
    //     BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //     std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //     for (int i = 0; i < size-1; ++i) {
    //         solutionVec[i] = solutionVec1[i];
    //     }
    //     auto solutionPrev = spacePrev->generateFunction(solutionVec);
    //     BEM::plotFunction("solutionOK"+std::to_string(order-50), *solutionOK, meshOK);
    //     BEM::plotFunction("solutionDerOK", solutionDerOKL, meshOK);
        
    //     for (int s = 1; s <= power-1; ++s) {
    //         size = std::pow(2, s);
    //         MeshCurve1D *mesh = new MeshCurve1D(size, curve);
    //         RegularP1_0Mesh_1D *space = new RegularP1_0Mesh_1D(*mesh);
    //         RiemannLiouvilleMesh rlm(*space, RiemannLiouvilleMesh::Side::LEFT, order, d, q);
    //         auto &matrix = rlm.getMatrix();
    //         Eigen::BDCSVD<BEM::Matrix, Eigen::ComputeThinU> SVD(matrix);
    //         auto vec = space->testAgainstBasis(rhs);
    //         BEM::ColVector solutionVec1(rlm.getMatrix().colPivHouseholderQr().solve(vec));
    //         std::vector<BEM::Complex> solutionVec(size-1, 0.0);
    //         for (int i = 0; i < size-1; ++i) {
    //             solutionVec[i] = solutionVec1[i];
    //         }
    //         auto solution = space->generateFunction(solutionVec);
    //         auto &solutionDerL = solution->derivative(order);
    //         // BEM::plotFunction("sol"+std::to_string(s), *solution, *mesh);
    //         // BEM::plotFunction("der"+std::to_string(s), solutionDerL, *mesh);
    //         std::cout << solution->L2Error(ExplicitScalarFunction_2D([&](double t, double s){return (*solutionOK)(t);})) << " " <<
    //                   solutionDerL.L2Error(ExplicitScalarFunction_2D([&](double t, double s){return solutionDerOKL(t);})) << std::endl;
    //         // delete mesh;
    //         // delete space;
    //         delete meshPrev;
    //         delete spacePrev;
    //         meshPrev = mesh;
    //         spacePrev = space;
    //     }
    // }

    
    // BOOST_AUTO_TEST_CASE(MatrixConstruction, *boost::unit_test::disabled())
    // {
    //     msg(1) << "MatrixConstruction" << endMsg;
    //     double error = 0;
    //     double prevError = 0;
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     for(int i = 0; i < 3; ++i) {            
    //         GeoP1_1D space(curve, std::pow(2, i)*10, 1./0.4, true);
    //         RiemannLiouville rl(space, RiemannLiouville::Side::LEFT, 1.4, d);
    //         rl.assembleMassMatrix();
    //         ExplicitScalarFunction_1D rhs([](double t)->BEM::Complex {return 1.0;});
    //         auto rhsVec = space.testAgainstBasis(rhs);
    //         BEM::ColVector solutionVec(rl.getMatrix().colPivHouseholderQr().solve(rhsVec));
        
    //         auto solution = space.DiscreteSpaceOnCurve_1D::generateFunction(solutionVec);
    //         ExplicitScalarFunction_1D diff([&solution](double t) {return (*solution)(t) - (-std::pow(t, 1.4) + std::pow(t, 0.4))/(1.4*std::tgamma(1.4));});
    //         prevError = error;
    //         error = BEM::L2Norm(diff);
    //         if (prevError > 0) {
    //             BOOST_CHECK(std::log(error/prevError)/std::log(0.5) > 1.7);
    //         }
    //     }
    // }

    // BOOST_AUTO_TEST_CASE(MatrixL2Construction, *boost::unit_test::disabled())
    // {
    //     msg(1) << "MatrixL2Construction" <<endMsg;
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     RegularP1_1D space(curve, 10, true);
    //     ExplicitScalarFunction_1D function([](double t) {return 11.0 + t*0;});
    //     L2 sb(function, space, space);
    //     sb.assembleMassMatrix();
    //     for (int i = 0; i < 10; ++i) {
    //         for (int j = 0; j < 10; ++j) {
    //             if (i == j) {
    //                 BOOST_CHECK_CLOSE(std::real(sb.getMatrix()(i, j)), 2./3., tolerance);
    //             } else if (std::abs(i - j) == 1){
    //                 BOOST_CHECK_CLOSE(std::real(sb.getMatrix()(i, j)), 1./6., tolerance);
    //             } else {
    //                 BOOST_CHECK(std::abs(sb.getMatrix()(i, j)) < tolerance);
    //             }
    //         }
    //     }
    // }
    BOOST_AUTO_TEST_SUITE_END()
    
    BOOST_AUTO_TEST_SUITE(FracParabolic)

    BOOST_AUTO_TEST_CASE(ConstructorTest)
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
        unsigned size = std::pow(2, 8);
        int order = 75;
        MeshCurve1D mesh(size, curve);
        RegularP1_L0Mesh_1D spaceTrial(mesh);
        RegularP0Mesh_1D spaceTest(mesh);
        FractionalParabolic rlm(spaceTrial, spaceTest, order, d);
    }

    BOOST_AUTO_TEST_CASE(MatrixTest)
    {
        TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
        ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
        unsigned size = std::pow(2, 3);
        int order = 90;
        MeshCurve1D mesh(size, curve);
        RegularP1_L0Mesh_1D spaceTrial(mesh);
        RegularP0Mesh_1D spaceTest(mesh);
        FractionalParabolic rlm(spaceTrial, spaceTest, order, d);
        rlm.assembleMassMatrix();
        ExplicitScalarFunction_2D f([](double t, double s)->BEM::Complex {return 1.0*t*t;});
        auto vec = spaceTest.testAgainstBasis(f);
        BEM::ColVector sol = rlm.getMatrix().colPivHouseholderQr().solve(vec);
        auto solFun = spaceTrial.generateFunction(BEM::toVector(sol));
        auto &dFun = solFun->derivative(order);
        ExplicitScalarFunction_2D f2([&](double t, double s)->BEM::Complex {return dFun(t);});
        auto vec2 = spaceTest.testAgainstBasis(f2);
        std::cout << vec2 << std::endl;
        BEM::plotFunction("testSol", *solFun, mesh);
        BEM::plotFunction("testSolDer", dFun, mesh);
    }

    BOOST_AUTO_TEST_SUITE_END()

    // BOOST_AUTO_TEST_SUITE(L2Time)

    // BOOST_AUTO_TEST_CASE(ConstructorTest)
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     unsigned size = std::pow(2, 8);
    //     MeshCurve1D mesh(size, curve);
    //     RegularP1_L0Mesh_1D spaceTrial(mesh);
    //     RegularP0Mesh_1D spaceTest(mesh);
    //     L2Mesh l2(spaceTrial, spaceTest, d);
    // }

    // BOOST_AUTO_TEST_CASE(MatrixTest)
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     unsigned size = std::pow(2, 2);
    //     MeshCurve1D mesh(size, curve);
    //     RegularP1_L0Mesh_1D spaceTrial(mesh);
    //     RegularP0Mesh_1D spaceTest(mesh);
    //     L2Mesh l2(spaceTrial, spaceTest, d);
    //     l2.assembleMassMatrix();
    //     ExplicitScalarFunction_2D f([](double t, double s)->BEM::Complex {return 1.0*t*t;});
    //     auto vec = spaceTest.testAgainstBasis(f);
    //     std::cout << vec << std::endl;
    //     BEM::ColVector sol = l2.getMatrix().colPivHouseholderQr().solve(vec);
    //     auto solFun = spaceTrial.generateFunction(BEM::toVector(sol));
    //     std::cout << sol << std::endl;
    //     BEM::plotFunction("testSol", *solFun, mesh);
    // }

    // BOOST_AUTO_TEST_CASE(SpaceTime)
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     TrigonometricCurve curveT(0.5, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return .0;});
    //     unsigned sizeTime = std::pow(2, 6);
    //     unsigned sizeSpace = std::pow(2, 5);
    //     MeshCurve1D meshTime(sizeTime, curveT);
    //     MeshCurve1D meshSpace(sizeSpace, curve);
    //     RegularP1_L0Mesh_1D spaceTrialTime(meshTime);
    //     RegularP0Mesh_1D spaceTestTime(meshTime);
    //     RegularP1_0Mesh_1D spaceSpace(meshSpace);
    //     int orderTime = BEM::getEnv<int>("FRACORDER");
    //     ExplicitScalarFunction_1D dRl([orderTime](double t)->BEM::Complex {return 1;});
    //     // Time
    //     FractionalParabolic frac(spaceTrialTime, spaceTestTime, orderTime, d);
    //     L2Mesh l2(spaceTrialTime, spaceTestTime, d);
    //     l2.assembleMassMatrix();
    //     frac.assembleMassMatrix();

    //     // Space
    //     int orderSpace = 99;
    //     RiemannLiouvilleMesh rlm(spaceSpace, RiemannLiouvilleMesh::Side::LEFT, orderSpace, dRl, q);
    //     L2Mesh l2s(spaceSpace, spaceSpace, d);

    //     rlm.assembleMassMatrix();
    //     l2s.assembleMassMatrix();
    //     std::unique_ptr<BEM::Matrix> _massMatrix(nullptr);
    //     _massMatrix.reset(new BEM::Matrix(spaceSpace.getSize()*spaceTrialTime.getSize(),spaceSpace.getSize()*spaceTrialTime.getSize()));
    //     std::unique_ptr<BEM::ColVector> _rhs(nullptr);
    //     _rhs.reset(new BEM::ColVector(spaceSpace.getSize()*spaceTrialTime.getSize()));
    //     for (int i = 0; i < spaceSpace.getSize(); ++i){
    //         for (int j = 0; j < spaceTestTime.getSize(); ++j){
    //             for (int k = 0; k < spaceSpace.getSize(); ++k){
    //                 for (int l = 0; l < spaceTestTime.getSize(); ++l){
    //                     _massMatrix->operator()(i + j*spaceSpace.getSize(), k + l*spaceSpace.getSize()) =
    //                         l2s.getMatrix()(i, k)*frac.getMatrix()(j, l) + rlm.getMatrix()(i, k)*l2.getMatrix()(j, l);
    //                 }
    //             }
    //         }            
    //     }
    //     ExplicitScalarFunction_2D fTime([](double t, double s)->BEM::Complex {
    //         if (t>0.25) {
    //             return 0.0;
    //         }
    //         return std::sin(M_PI*t/0.25);
    //     });
    //     auto vecTime = spaceTestTime.testAgainstBasis(fTime);
        
    //     ExplicitScalarFunction_2D fSpace([](double t, double s)->BEM::Complex {return t*(1-t);});
    //     auto vecSpace = spaceSpace.testAgainstBasis(fSpace);
        
    //     for (int i = 0; i < spaceSpace.getSize(); ++i){
    //         for (int j = 0; j < spaceTestTime.getSize(); ++j){
    //             _rhs->operator()(i+j*spaceSpace.getSize())=vecTime[j]*vecSpace[i];
    //         }
    //     }

    //     BEM::ColVector fullSol = _massMatrix->colPivHouseholderQr().solve(*_rhs);
    //     auto vecSol = BEM::toVector(fullSol);

    //     for (int j = 0 ; j < spaceTrialTime.getSize(); ++j) {
    //         std::vector<BEM::Complex> vec(vecSol.begin()+j*spaceSpace.getSize(), vecSol.begin()+(j+1)*spaceSpace.getSize());
    //         auto solFun = spaceSpace.generateFunction(vec);
    //         BEM::plotFunction("testSol"+std::to_string(orderTime)+"_"+std::to_string(j), *solFun, meshSpace);
    //     }
    // }

    // BOOST_AUTO_TEST_CASE(SpaceTimePara)
    // {
    //     TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     TrigonometricCurve curveT(1, 0, std::vector<double>{0}, std::vector<double>{0});
    //     ExplicitScalarFunction_1D d([](double t)->BEM::Complex {return 1.0;});
    //     ExplicitScalarFunction_1D q([](double t)->BEM::Complex {return .0;});
    //     unsigned sizeTime = std::pow(2, 6);
    //     unsigned sizeSpace = std::pow(2, 5);
    //     MeshCurve1D meshTime(sizeTime, curveT);
    //     MeshCurve1D meshSpace(sizeSpace, curve);
    //     RegularP1_L0Mesh_1D spaceTrialTime(meshTime);
    //     RegularP0Mesh_1D spaceTestTime(meshTime);
    //     RegularP1_0Mesh_1D spaceSpace(meshSpace);
    //     int orderTime = BEM::getEnv<int>("FRACORDER");
    //     // Time
    //     FractionalParabolic frac(spaceTrialTime, spaceTestTime, orderTime, d);
    //     L2Mesh l2(spaceTrialTime, spaceTestTime, d);
    //     l2.assembleMassMatrix();
    //     frac.assembleMassMatrix();
    //     std::unique_ptr<BEM::Matrix> svdMat(nullptr);
    //     svdMat.reset(new BEM::Matrix(spaceSpace.getSize(), spaceTrialTime.getSize()*20));
    //     int colcounter = 0;
    //     std::uniform_real_distribution<double> _unif(0.1, 1.9);
    //     std::mt19937 _rng((std::random_device())());


    //     for (int rand = 0; rand < 10; rand++) {
    //         // Space
    //         int orderSpace = 90;
    //         double a = _unif(_rng);
    //         double b = _unif(_rng);
    //         double c = _unif(_rng);
    //         ExplicitScalarFunction_1D dRl([&](double t)->BEM::Complex {
    //             if (t < 0.25) {
    //                 return 1.;
    //             }
    //             if (t < 0.5) {
    //                 return a;
    //             }
    //             if (t < 0.75) {
    //                 return b;
    //             }
    //             return c;
    //         });
    //         RiemannLiouvilleMesh rlm(spaceSpace, RiemannLiouvilleMesh::Side::LEFT, orderSpace, dRl, q);
    //         L2Mesh l2s(spaceSpace, spaceSpace, d);

    //         rlm.assembleMassMatrix();
    //         l2s.assembleMassMatrix();
    //         std::unique_ptr<BEM::Matrix> _massMatrix(nullptr);
    //         _massMatrix.reset(new BEM::Matrix(spaceSpace.getSize()*spaceTrialTime.getSize(),spaceSpace.getSize()*spaceTrialTime.getSize()));
    //         std::unique_ptr<BEM::ColVector> _rhs(nullptr);
    //         _rhs.reset(new BEM::ColVector(spaceSpace.getSize()*spaceTrialTime.getSize()));
    //         for (int i = 0; i < spaceSpace.getSize(); ++i){
    //             for (int j = 0; j < spaceTestTime.getSize(); ++j){
    //                 for (int k = 0; k < spaceSpace.getSize(); ++k){
    //                     for (int l = 0; l < spaceTestTime.getSize(); ++l){
    //                         _massMatrix->operator()(i + j*spaceSpace.getSize(), k + l*spaceSpace.getSize()) =
    //                             l2s.getMatrix()(i, k)*frac.getMatrix()(j, l) + rlm.getMatrix()(i, k)*l2.getMatrix()(j, l);
    //                     }
    //                 }
    //             }            
    //         }
    //         ExplicitScalarFunction_2D fTime([](double t, double s)->BEM::Complex {
    //             return std::sin(M_PI*(t)/0.1)+1.;
    //         });
    //         auto vecTime = spaceTestTime.testAgainstBasis(fTime);
        
    //         ExplicitScalarFunction_2D fSpace([](double t, double s)->BEM::Complex {return t*(1-t);});
    //         auto vecSpace = spaceSpace.testAgainstBasis(fSpace);
        
    //         for (int i = 0; i < spaceSpace.getSize(); ++i){
    //             for (int j = 0; j < spaceTestTime.getSize(); ++j){
    //                 _rhs->operator()(i+j*spaceSpace.getSize())=vecTime[j]*vecSpace[i];
    //             }
    //         }

    //         BEM::ColVector fullSol = _massMatrix->colPivHouseholderQr().solve(*_rhs);
    //         auto vecSol = BEM::toVector(fullSol);

    //         for (int j = 0 ; j < spaceTrialTime.getSize(); ++j) {
    //             std::vector<BEM::Complex> vec(vecSol.begin()+j*spaceSpace.getSize(), vecSol.begin()+(j+1)*spaceSpace.getSize());
    //             for (int k = 0; k < spaceSpace.getSize(); ++k) {
    //                 svdMat->operator()(k,colcounter) = vec[k];
    //             }
    //             auto solFun = spaceSpace.generateFunction(vec);
    //             BEM::plotFunction("testSol"+std::to_string(orderTime)+"_"+std::to_string(j), *solFun, meshSpace);
    //             colcounter +=1;
    //         }
    //         std::cout << "done with k = " << rand << std::endl;
    //     }
    //     Eigen::BDCSVD<BEM::Matrix, Eigen::ComputeThinU> SVD(*svdMat);
    //     auto _singVals = SVD.singularValues();
    //     std::cout << _singVals << std::endl;
    // }
    
    // BOOST_AUTO_TEST_SUITE_END()

}

#endif
