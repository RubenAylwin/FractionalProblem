#include <Curve.h>
#include <Point.h>
#include <ScalarValuedFunction.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

/**
 * @brief: Get point at t (depends on given parametrization).
 * @desc: Wrapper to specific implementation.
 */
Point2D Curve2D::at(const double t) const
{
    assert(t >= _lowLimit - 0.00001 and t <= _uppLimit + 0.00001);
    return evaluateAt(t);
}

/**
 * @brief: Get normal at point t.
 */
Vector2D Curve2D::normal(const double t) const
{
    assert(t >= _lowLimit and t <= _uppLimit);
    return evaluateNormalAt(t);
}

/**
 * @brief: Get jacobian at point t.
 */
double Curve2D::jacobian(const double t) const
{
    return normal(t).norm();
}

/**
 * @brief: Get jacobian as a function.
 */
std::unique_ptr<ScalarFunctionBase_1D> Curve2D::surfaceMeasure(void) const
{
    return std::unique_ptr<ScalarFunctionBase_1D>(new ExplicitScalarFunction_1D([&](double t){return this->jacobian(t);}));
}

/**
 * @brief: Constructor. Depends only on desired length.
 */
StraightCurve::StraightCurve(double length) :
    Curve2D(0, length),
    _length{length}
{
}

/**
 * @brief: Specific implementation of evaluation.
 */
Point2D StraightCurve::evaluateAt(const double t) const
{
    return Point2D(t, 0.);
}

/**
 * @brief: Specific implementation of normal..
 */
Vector2D StraightCurve::evaluateNormalAt(const double t) const
{
    return Vector2D(t*0., 1.0);
}

/**
 * @brief: Constructor. Requires period, heigth and sine & cosine coefficients.
 */
TrigonometricCurve::TrigonometricCurve(double period, double height, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients) :
    PeriodicCurve(period),
    _height{height},
    _sineCoefficients{sineCoefficients},
    _cosineCoefficients{cosineCoefficients}
{
}

/**
 * @brief: Returns the parameters (copy).
 */
std::vector<double> TrigonometricCurve::getParameters(void) const
{
    std::vector<double> parameters;
    parameters.insert(parameters.end(), _sineCoefficients.begin(), _sineCoefficients.end());
    parameters.insert(parameters.end(), _cosineCoefficients.begin(), _cosineCoefficients.end());
    return parameters;
}

/**
 * @brief: Specific implementation of evaluation.
 */
Point2D TrigonometricCurve::evaluateAt(const double t) const
{
    double X = _period*t;
    double Y = _height;
    for (size_t i = 1 ; i <= _sineCoefficients.size(); ++i) {
        Y += _sineCoefficients[i - 1]*std::sin(2.*M_PI*i*t);
    }
    for (size_t i = 1 ; i <= _cosineCoefficients.size(); ++i) {
        Y += _cosineCoefficients[i - 1]*std::cos(2.*M_PI*i*t);
    }
    return Point2D(X, Y);
}

/**
 * @brief: Specific implementation of normal.
 */
Vector2D TrigonometricCurve::evaluateNormalAt(const double t) const
{
    double Y = _period;
    double X = 0;
    for (size_t i = 1 ; i <= _sineCoefficients.size(); ++i) {
        X += -1.*_sineCoefficients[i - 1]*std::cos(2.*M_PI*i*t)*2.*M_PI*i;
    }
    for (size_t i = 1 ; i <= _cosineCoefficients.size(); ++i) {
        X += _cosineCoefficients[i - 1]*std::sin(2.*M_PI*i*t)*2.*M_PI*i;
    }

    return Vector2D(X, Y);
}

/**
 * @brief: Constructor. Built from an ordered list of points (vertices).
 */
PolyPeriodicCurve::PolyPeriodicCurve(double period, std::vector<Point2D> &&orderedPoints) :
    PeriodicCurve(period),
    _numLines(orderedPoints.size() - 1),
    _orderedPoints(orderedPoints)
{
    //First and last point should "meet"
    assert(_orderedPoints.front().getX() == 0);
    assert(_orderedPoints.back().getX() == period);
    assert(_orderedPoints.front().getY() == _orderedPoints.back().getY());
    for (unsigned i = 0; i <= _numLines; ++i) {
        _intervalPartition.push_back(i*1.0/_numLines);
    }

}

/**
 * @brief: Specific implementation of evaluation.
 */
Point2D PolyPeriodicCurve::evaluateAt(const double t) const
{
    if (t == 0) {
        return _orderedPoints.front();
    }

    if (t == 1) {
        return _orderedPoints.back();
    }
    auto upperBound = std::upper_bound(_intervalPartition.begin(), _intervalPartition.end(), t);
    auto upperIndex = upperBound - _intervalPartition.begin();
    int lowerIndex = upperIndex - 1;
    return Point2D::interpolate(_orderedPoints[lowerIndex], _orderedPoints[upperIndex], (t - _intervalPartition[lowerIndex])/(_intervalPartition[upperIndex] - _intervalPartition[lowerIndex]));
    
}

/**
 * @brief: Specific implementation of normal.
 */
Vector2D PolyPeriodicCurve::evaluateNormalAt(const double t) const
{
    auto upperBound = std::upper_bound(_intervalPartition.begin(), _intervalPartition.end(), t);
    auto upperIndex = upperBound - _intervalPartition.begin();
    int lowerIndex = upperIndex - 1;
    auto oPoint = _orderedPoints[lowerIndex];
    auto dPoint = _orderedPoints[upperIndex];
    auto diffPoint = dPoint - oPoint;
    double Jac = (_intervalPartition[upperIndex] - _intervalPartition[lowerIndex]);
    double Y = diffPoint.getX()/(Jac);
    double X = -diffPoint.getY()/(Jac);
    return Vector2D(X, Y);
}
