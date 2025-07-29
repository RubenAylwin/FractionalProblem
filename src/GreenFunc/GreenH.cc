#include <GreenH.h>
#include <Utilities.h>
#include <Point.h>
#include <iostream>
#include <cassert>
#include <cmath>

/**
 * @brief: Constructor
 */
GreenH2D::GreenH2D(double wavenumber) :
    _wavenumber(wavenumber)
{
    assert(wavenumber > 0 and "Error: wavenumber should be positive");
}

/**
 * @brief: Evaluation on 2D points
 */
BEM::Complex GreenH2D::operator()(const Point2D &X, const Point2D &Y) const
{
    return BEM::hankel0_1(_wavenumber*(X - Y).norm())*(-BEM::I/4.0);
}
