#include <GreenL.h>
#include <cassert>
#include <cmath>
#include <Point.h>
#include <iostream>

/**
 * @brief: Evaluation on 2D points.
 */
BEM::Complex GreenL2D::operator()(const Point2D &X, const Point2D &Y) const
{
    return std::log((X - Y).norm())*1./(2.*M_PI);
}


