#ifndef MESH_TEST
#define MESH_TEST
/*
 *TESTS FOR GEOMETRY CLASSES AND FUNCTIONS
 */

#include <Point.h>
#include <Curve.h>
#include <Mesh.h>

namespace meshTest {
static double tolerance = 1E-9;
BOOST_AUTO_TEST_SUITE(Mesh)
BOOST_AUTO_TEST_SUITE(Straight)

BOOST_AUTO_TEST_CASE(Constructor)
{
    StraightCurve curve(10);
    MeshCurve1D mesh(20, curve);
}

BOOST_AUTO_TEST_CASE(Elements)
{
    StraightCurve curve(10);
    MeshCurve1D mesh(17, curve);
    auto el0 = mesh.getElement(0);
    BOOST_CHECK_CLOSE(el0.getA().getX(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getA().getY(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getX(), 10./17., tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getY(), 0.0, tolerance);

    auto el5 = mesh.getElement(5);
    BOOST_CHECK_CLOSE(el5.getA().getX(), 50./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getA().getY(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getX(), 60./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getY(), 0.0, tolerance);

}

BOOST_AUTO_TEST_CASE(Point)
{
    StraightCurve curve(10);
    MeshCurve1D mesh(17, curve);
    auto point = mesh.point(5);
    BOOST_CHECK_CLOSE(point.getX(), 50./17., tolerance);
    auto point2 = mesh.point(7);
    BOOST_CHECK_CLOSE(point2.getX(), 70./17., tolerance);
}
 
BOOST_AUTO_TEST_CASE(ElementFromPoint)
{
    StraightCurve curve(10);
    MeshCurve1D mesh(17, curve);
    auto el0 = mesh.elementFromBPoint(1);
    BOOST_CHECK_CLOSE(el0.getA().getX(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getA().getY(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getX(), 10./17., tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getY(), 0.0, tolerance);

    auto el5 = mesh.elementFromAPoint(5);
    BOOST_CHECK_CLOSE(el5.getA().getX(), 50./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getA().getY(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getX(), 60./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getY(), 0.0, tolerance);
}
 
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Trigonometric)

BOOST_AUTO_TEST_CASE(Constructor)
{
    TrigonometricCurve curve(10, 0, std::vector<double>{1.0}, std::vector<double>{0.0});
    MeshCurve1D mesh(20, curve);
}

BOOST_AUTO_TEST_CASE(Elements)
{
    TrigonometricCurve curve(10, 0, std::vector<double>{1.0}, std::vector<double>{0.0});
    MeshCurve1D mesh(17, curve);
    auto el0 = mesh.getElement(0);
    BOOST_CHECK_CLOSE(el0.getA().getX(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getA().getY(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getX(), 10./17., tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getY(), std::sin(2.0*M_PI*10./(17.*10.)), tolerance);

    auto el5 = mesh.getElement(5);
    BOOST_CHECK_CLOSE(el5.getA().getX(), 50./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getA().getY(), std::sin(2.0*M_PI*50./(17.*10.)), tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getX(), 60./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getY(), std::sin(2.0*M_PI*60./(17.*10.)), tolerance);

}

BOOST_AUTO_TEST_CASE(Point)
{
    TrigonometricCurve curve(10, 0, std::vector<double>{1.0}, std::vector<double>{0.0});
    MeshCurve1D mesh(17, curve);
    auto point = mesh.point(5);
    BOOST_CHECK_CLOSE(point.getX(), 50./17., tolerance);
    BOOST_CHECK_CLOSE(point.getY(), std::sin(2.0*M_PI*50./(17.*10.)), tolerance);
    auto point2 = mesh.point(7);
    BOOST_CHECK_CLOSE(point2.getX(), 70./17., tolerance);
    BOOST_CHECK_CLOSE(point2.getY(), std::sin(2.0*M_PI*70./(17.*10.)), tolerance);
}
 
BOOST_AUTO_TEST_CASE(ElementFromPoint)
{
    TrigonometricCurve curve(10, 0, std::vector<double>{1.0}, std::vector<double>{0.0});
    MeshCurve1D mesh(17, curve);
    auto el0 = mesh.elementFromBPoint(1);
    BOOST_CHECK_CLOSE(el0.getA().getX(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getA().getY(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getX(), 10./17., tolerance);
    BOOST_CHECK_CLOSE(el0.getB().getY(), std::sin(2.0*M_PI*10./(17.*10.)), tolerance);
    auto el5 = mesh.elementFromAPoint(5);
    BOOST_CHECK_CLOSE(el5.getA().getX(), 50./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getA().getY(), std::sin(2.0*M_PI*50./(17.*10.)), tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getX(), 60./17., tolerance);
    BOOST_CHECK_CLOSE(el5.getB().getY(), std::sin(2.0*M_PI*60./(17.*10.)), tolerance);}
 
BOOST_AUTO_TEST_CASE(ElementContainingPoint)
{
    TrigonometricCurve curve(1, 0, std::vector<double>{1.0}, std::vector<double>{0.0});
    MeshCurve1D mesh(17, curve);
    BOOST_CHECK_EQUAL(mesh.elementWithPoint(0.).first, 0);
    BOOST_CHECK_EQUAL(mesh.elementWithPoint(1./17.).first, 0);
    BOOST_CHECK_EQUAL(mesh.elementWithPoint(1./17. + tolerance).first, 1);
    BOOST_CHECK_EQUAL(mesh.elementWithPoint(0.28).first, 4);
    BOOST_CHECK_EQUAL(mesh.elementWithPoint(0.38).first, 6);
    BOOST_CHECK_EQUAL(mesh.elementWithPoint(1.).first, 16);
}
    
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
    }     
#endif
