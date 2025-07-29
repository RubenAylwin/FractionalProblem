#include <GeometricVector.h>
#include <cmath>
#include <iostream>


/**
 * @brief: Constructor.
 */
Vector2D::Vector2D(double x, double y) :
    _x{x},
    _y{y}
{ }

/**
 * @brief: 2-Norm.
 */
double Vector2D::norm(void) const {
    return std::sqrt(_x*_x + _y*_y);
}

/**
 * @brief: To print the vector.
 * TODO: overload <<
 */
void Vector2D::show(void) const {
    std::cout << "Vector[" << _x << "," << _y << "]" << std::endl;
}
