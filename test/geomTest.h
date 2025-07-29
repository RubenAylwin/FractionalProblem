#ifndef GEOM_TEST
#define GEOM_TEST
/*
 *TESTS FOR GEOMETRY CLASSES AND FUNCTIONS
 */
#include <Point.h>
#include <GeometricVector.h>
#include <Curve.h>
#include <cmath>
#include <iostream>
static double tolerance = 1E-9;

BOOST_AUTO_TEST_SUITE(Geometry)

/*
 *TEST FOR POINT
 */
BOOST_AUTO_TEST_SUITE(POINT)

BOOST_AUTO_TEST_CASE(ConstructorTest)
{
    msg(1) << "start Geometry::Point::ConstructorTest" << endMsg;    
    Point2D p1(0, 1);
    BOOST_CHECK(p1.getX() == 0);
    BOOST_CHECK(p1.getY() == 1);
    msg(1) << "end Geometry::Point::ConstructorTest" << endMsg;    
}

BOOST_AUTO_TEST_CASE(AdditionTest)
{
    msg(1) << "start Geometry::Point::AdditionTest" << endMsg;    
    Point2D p1(0, 1);
    Point2D p2(1, 0);
    Point2D p3 = p1 + p2;
    BOOST_CHECK(p3.getX() == 1);
    BOOST_CHECK(p3.getY() == 1);
    msg(1) << "end Geometry::Point::AdditionTest" << endMsg;    
}

BOOST_AUTO_TEST_CASE(SubtractionTest)
{
    msg(1) << "start Geometry::Point::SubtractionTest" << endMsg;    
    Point2D p1(0, 1);
    Point2D p2(1, 0);
    Point2D p3 = p1 - p2;
    BOOST_CHECK(p3.getX() == -1);
    BOOST_CHECK(p3.getY() == 1);
    msg(1) << "end Geometry::Point::SubtractionTest" << endMsg;    
}
     
BOOST_AUTO_TEST_SUITE_END()
/*
 *END TEST FOR POINT
 */

/*
 *TEST FOR VECTOR
 */
BOOST_AUTO_TEST_SUITE(VECTOR)

BOOST_AUTO_TEST_CASE(ConstructorTest)
{
    msg(1) << "start Geometry::Vector::ConstructorTest" << endMsg;    
    Vector2D v1(0, 1);
    BOOST_CHECK(v1.getX() == 0);
    BOOST_CHECK(v1.getY() == 1);
    msg(1) << "end Geometry::Vector::ConstructorTest" << endMsg;    
}

BOOST_AUTO_TEST_CASE(AdditionTest)
{
    msg(1) << "start Geometry::Vector::AdditionTest" << endMsg;    
    Vector2D v1(0, 1);
    Vector2D v2(1, 0);
    Vector2D v3 = v1 + v2;
    BOOST_CHECK(v3.getX() == 1);
    BOOST_CHECK(v3.getY() == 1);
    msg(1) << "end Geometry::Vector::AdditionTest" << endMsg;    
}

BOOST_AUTO_TEST_CASE(SubtractionTest)
{
    msg(1) << "start Geometry::Vector::SubtractionTest" << endMsg;    
    Vector2D v1(0, 1);
    Vector2D v2(1, 0);
    Vector2D v3 = v1 - v2;
    BOOST_CHECK(v3.getX() == -1);
    BOOST_CHECK(v3.getY() == 1);
    msg(1) << "end Geometry::Vector::SubtractionTest" << endMsg;    
}

BOOST_AUTO_TEST_CASE(ScalarProductTest)
{
    msg(1) << "start Geometry::Vector::ScalarProductTest" << endMsg;    
    Vector2D v1(0, 1);
    Vector2D v2 = v1*3;
    BOOST_CHECK(v2.getX() == 0);
    BOOST_CHECK(v2.getY() == 3);
    msg(1) << "end Geometry::Vector::ScalarProductTest" << endMsg;    
}


BOOST_AUTO_TEST_SUITE_END()
/*
 *END TEST FOR VECTOR
 */

/*
 *TEST FOR TRIGONOMETRIC CURVE
 */
BOOST_AUTO_TEST_SUITE(TrigonometricCurveTest)

BOOST_AUTO_TEST_CASE(ConstructAndEvaluate)
{
    msg(1) << "start Geometry::TrigonometricCurveTest::ConstructAndEvaluate " << endMsg;
    TrigonometricCurve c1(1, 0, std::vector<double>{0}, std::vector<double>{0});
    auto p = c1.at(0.2);
    auto n = c1.normal(0.3);
    BOOST_CHECK_CLOSE(p.getX(), 0.2, tolerance);
    BOOST_CHECK_CLOSE(p.getY(), 0, tolerance);
    BOOST_CHECK_CLOSE(n.getX(), 0, tolerance);
    msg(1) << "end Geometry::TrigonometricCurveTest::ConstructAndEvaluate " << endMsg;
}

BOOST_AUTO_TEST_CASE(ConstructAndEvaluate2)
{
    msg(1) << "start Geometry::TrigonometricCurveTest::ConstructAndEvaluate2 " << endMsg;
    TrigonometricCurve c1(1, 0, std::vector<double>{3}, std::vector<double>{0});
    auto p = c1.at(0.2);
    auto n = c1.normal(0.3);

    double nx = -3*2*M_PI*std::cos(2*M_PI*0.3);
    double ny = 1;
    double norm = std::sqrt(nx*nx+1);
    
    BOOST_CHECK_CLOSE(p.getX(), 0.2, tolerance);
    BOOST_CHECK_CLOSE(p.getY(), std::sin(2*M_PI*0.2)*3, tolerance);
    BOOST_CHECK_CLOSE(n.norm(), norm, tolerance);
    BOOST_CHECK_CLOSE(n.getX(), nx, tolerance);
    BOOST_CHECK_CLOSE(n.getY(), ny, tolerance);
    msg(1) << "end Geometry::TrigonometricCurveTest::ConstructAndEvaluate2 " << endMsg;
}

BOOST_AUTO_TEST_CASE(ConstructAndEvaluate3)
{
    msg(1) << "start Geometry::TrigonometricCurveTest::ConstructAndEvaluate3 " << endMsg;
    TrigonometricCurve c1(1, 0, std::vector<double>{3, 0, 0, 2}, std::vector<double>{0, 5, 1});
    auto n = c1.normal(0.3);
    double nx = -3*2*M_PI*std::cos(2*M_PI*0.3) - 2*2*M_PI*4*std::cos(2*M_PI*0.3*4) + 5*2*M_PI*2*std::sin(2*M_PI*0.3*2) + 2*M_PI*3*std::sin(2*M_PI*0.3*3);
    double ny = 1;
    double norm = std::sqrt(nx*nx+1);
    BOOST_CHECK_CLOSE(n.norm(), norm, tolerance);
    BOOST_CHECK_CLOSE(n.getX(), nx, tolerance);
    BOOST_CHECK_CLOSE(n.getY(), ny, tolerance);
    msg(1) << "end Geometry::TrigonometricCurveTest::ConstructAndEvaluate3 " << endMsg;
}


BOOST_AUTO_TEST_CASE(JacobianTest) {
    msg(1) << "start Geometry::TrigonometricCurveTest::JacobianTest" << endMsg;
    TrigonometricCurve c(1, 0, std::vector<double>{3}, std::vector<double>{0});
    auto jacobian = c.surfaceMeasure();
    msg(1) << "end Geometry::TrigonometricCurveTest::JacobianTest" << endMsg;
}

     
BOOST_AUTO_TEST_SUITE_END()
/*
 *END TEST FOR TRIGONOMETRIC CURVE
 */

BOOST_AUTO_TEST_SUITE(PolyPerCurve)

BOOST_AUTO_TEST_CASE(ConstructAndEvaluateTriangle) {
    msg(1) << "start Geometry::PolyPerCurve::ConstructAndEvaluateTriangle " << endMsg;
    PolyPeriodicCurve curve(1, std::vector<Point2D>{Point2D(0,0), Point2D(0.5, 1), Point2D(1, 0)});
    {
        auto p = curve.at(0);
        BOOST_CHECK_CLOSE(p.getX(), 0, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 0, tolerance);
    }
    {
        auto p = curve.at(0.5);
        BOOST_CHECK_CLOSE(p.getX(), 0.5, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 1, tolerance);
    }
    {
        auto p = curve.at(1);
        BOOST_CHECK_CLOSE(p.getX(), 1, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 0, tolerance);
    }
    {
        auto p = curve.at(.25);
        BOOST_CHECK_CLOSE(p.getX(), .25, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), .5, tolerance);
    }
    {
        auto p = curve.at(.75);
        BOOST_CHECK_CLOSE(p.getX(), .75, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), .5, tolerance);
    }
    msg(1) << "end Geometry::PolyPerCurve::ConstructAndEvaluateTriangle " << endMsg;    
}

BOOST_AUTO_TEST_CASE(ConstructAndEvaluatePolyhedra) {
    PolyPeriodicCurve curve(1, std::vector<Point2D>{Point2D(0,0), Point2D(0.25, 0), Point2D(.25, .5), Point2D(.75, .5), Point2D(.75, 0), Point2D(1, 0)});
    msg(1) << "start Geometry::PolyPerCurve::ConstructAndEvaluatePolyhedra " << endMsg;
    {
        auto p = curve.at(0);
        BOOST_CHECK_CLOSE(p.getX(), 0, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 0, tolerance);
    }
    {
        auto p = curve.at(0.1);
        BOOST_CHECK_CLOSE(p.getX(), 0.125, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 0, tolerance);
    }
    {
        auto p = curve.at(.5);
        BOOST_CHECK_CLOSE(p.getX(), .5, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), .5, tolerance);
    }
    {
        auto p = curve.at(.6);
        BOOST_CHECK_CLOSE(p.getX(), .75, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), .5, tolerance);
    }

    {
        auto p = curve.at(.9);
        BOOST_CHECK_CLOSE(p.getX(), .875, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 0, tolerance);
    }

    {
        auto p = curve.at(.75);
        BOOST_CHECK_CLOSE(p.getX(), .75, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), .5/4, tolerance);
    }
    msg(1) << "end Geometry::PolyPerCurve::ConstructAndEvaluatePolyhedra " << endMsg;
}

BOOST_AUTO_TEST_CASE(NormalAtTriangle) {
    msg(1) << "start Geometry::PolyPerCurve::NormalAtTriangle " << endMsg;
    PolyPeriodicCurve curve(1, std::vector<Point2D>{Point2D(0,0), Point2D(0.5, 1), Point2D(1, 0)});
    {
        auto n = curve.normal(0);
        BOOST_CHECK_CLOSE(n.getX(), -2, tolerance);
        BOOST_CHECK_CLOSE(n.getY(), 1, tolerance);
    }
    {
        auto n = curve.normal(0.5);
        BOOST_CHECK_CLOSE(n.getX(), 2, tolerance);
        BOOST_CHECK_CLOSE(n.getY(), 1, tolerance);
    }
    {
        auto p = curve.normal(.9);
        BOOST_CHECK_CLOSE(p.getX(), 2, tolerance);
        BOOST_CHECK_CLOSE(p.getY(), 1, tolerance);
    }
    msg(1) << "end Geometry::PolyPerCurve::NormalAtTriangle " << endMsg;
}


BOOST_AUTO_TEST_SUITE_END()
/*
 *END TEST FOR POLY PER CURVE
 */


BOOST_AUTO_TEST_SUITE_END()
#endif
