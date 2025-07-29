#include <Point.h>
#include <cmath>
#include <iostream>
#include <cassert>

/**
 * @brief: Constructor.
 */
Point2D::Point2D(double x, double y) :
    _x{x},
    _y{y}
{
}

/**
 * @brief: Norm.
 */
double Point2D::norm(void) const {
    return std::sqrt(_x*_x + _y*_y);
}

/**
 * @brief: Print.
 */
void Point2D::show(void) const {
    std::cout << "Point[" << _x << "," << _y << "]" << std::endl;
}

/**
 * @brief: Interpolate between 2 points
 */
Point2D Point2D::interpolate(const Point2D &first, const Point2D &second, double t)
{
    assert(0 <= t);
    assert(t <= 1);
    return first*(1 - t) + second*t;
}
