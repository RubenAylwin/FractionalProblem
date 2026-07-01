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
        
        unsigned size = std::pow(2, 3);
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
            return std::abs(std::cos(M_PI*fOrder));
        };
        std::vector<BEM::Interval1D> limits{BEM::Interval1D(0., 5.), BEM::Interval1D(1., 1.),  BEM::Interval1D(1., 1.)};
        factory.trainGreedy(limits, 100, 1e-6, est);

    }
    
    BOOST_AUTO_TEST_SUITE_END()
    
}

#endif
