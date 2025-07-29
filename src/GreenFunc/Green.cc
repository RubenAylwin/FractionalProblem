#include <Green.h>
#include <Point.h>

/**
 * @brief: Smooth part of the green function.
 */
BEM::Complex GreenFunction2D::smoothPart(const Point2D &X, const Point2D &Y) const
{
    return operator()(X, Y) - singularity(X, Y);
}
