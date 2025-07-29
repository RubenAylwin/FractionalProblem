#include <TwoDimensionalIntegration.h>
#include <ScalarValuedFunction.h>
#include <stdexcept>


/**
 * @brief: Constructor. Takes number of points and limits for the construction of the points.
 */
Integrator_2D::Integrator_2D(unsigned n, double firstLowerLimitQuad, double firstUpperLimitQuad, double secondLowerLimitQuad, double secondUpperLimitQuad):
    _weights(n, 0),
    _points(n, Point2D(0, 0)),
    _firstLowerLimitQuad{firstLowerLimitQuad},
    _firstUpperLimitQuad(firstUpperLimitQuad),
    _secondLowerLimitQuad(secondLowerLimitQuad),
    _secondUpperLimitQuad(secondUpperLimitQuad)
{
}

/**
 * @brief: Integration of two dimensional function.
 */
BEM::Complex Integrator_2D::integrate(const double firstLowerLimit, const double firstUpperLimit, const double secondLowerLimit, const double secondUpperLimit, const ScalarFunctionBase_2D &integrand) const
{
    BEM::Complex result(0, 0);
    const double firstJacobian = ((firstUpperLimit - firstLowerLimit)/(_firstUpperLimitQuad - _firstLowerLimitQuad));
    const double secondJacobian = ((secondUpperLimit - secondLowerLimit)/(_secondUpperLimitQuad - _secondLowerLimitQuad));
    if (firstLowerLimit == secondLowerLimit and firstUpperLimit == secondUpperLimit) {
        for (unsigned i = 0; i < _weights.size(); ++i) {
            double evaluationPointX = firstLowerLimit + firstJacobian*(_points[i].getX() - _firstLowerLimitQuad);
            double secondJacobianFirst = secondJacobian*(1.0 - (_points[i].getX() - _firstLowerLimitQuad)/(_firstUpperLimitQuad-_firstLowerLimitQuad));
            double evaluationPointYFirst = secondLowerLimit + secondJacobianFirst*(_points[i].getY() - _secondLowerLimitQuad);

            double secondJacobianSecond = secondJacobian*((_points[i].getX() - _firstLowerLimitQuad)/(_firstUpperLimitQuad-_firstLowerLimitQuad));
            double evaluationPointYSecond = secondUpperLimit - (evaluationPointX - firstLowerLimit) + secondJacobianSecond*(_points[i].getY() - _secondLowerLimitQuad);
            result += _weights[i]*(integrand(evaluationPointX, evaluationPointYFirst)*secondJacobianFirst + integrand(evaluationPointX, evaluationPointYSecond)*secondJacobianSecond);
        }
        return result*firstJacobian;
    }

    for (unsigned i = 0; i < _weights.size(); ++i) {
        double evaluationPointX = firstLowerLimit + firstJacobian*(_points[i].getX() - _firstLowerLimitQuad);
        double evaluationPointY = secondLowerLimit + secondJacobian*(_points[i].getY() - _secondLowerLimitQuad);
        result += _weights[i]*integrand(evaluationPointX, evaluationPointY);
    }
    return firstJacobian*secondJacobian*result;
}

/**
 * @brief: Override for pointer function.
 */
BEM::Complex Integrator_2D::integrate(const double firstLowerLimit, const double firstUpperLimit, const double secondLowerLimit, const double secondUpperLimit, const std::unique_ptr<ScalarFunctionBase_2D> &integrand) const
{
    assert(integrand);
    return integrate(firstLowerLimit, firstUpperLimit, secondLowerLimit, secondUpperLimit, *integrand);
}

/**
 * @brief: Override for interval given as BEM::Interval1D.
 */
BEM::Complex Integrator_2D::integrate(const BEM::Interval1D firstInterval, const BEM::Interval1D secondInterval, const ScalarFunctionBase_2D &integrand) const
{
    return integrate(firstInterval.first, firstInterval.second, secondInterval.first, secondInterval.second, integrand);
}

/**
 * @brief: Override for interval given as BEM::Interval1D and pointer function.
 */
BEM::Complex Integrator_2D::integrate(const BEM::Interval1D firstInterval, const BEM::Interval1D secondInterval, const std::unique_ptr<ScalarFunctionBase_2D> &integrand) const
{
    return integrate(firstInterval, secondInterval, *integrand);
}

/**
 * @brief: Override for interval given as support list.
 */
BEM::Complex Integrator_2D::integrate(const BEM::Support1DL firstSupport, const BEM::Support1DL secondSupport, const ScalarFunctionBase_2D &integrand) const
{
    BEM::Complex result = .0;
    for (const auto &interval1 : firstSupport) {
        for (const auto &interval2 : secondSupport) {
            result += integrate(interval1, interval2, integrand);
        }
    }
    return result;
}

/**
 * @brief: Override for interval given as support list and pointer function.
 */
BEM::Complex Integrator_2D::integrate(const BEM::Support1DL firstSupport, const BEM::Support1DL secondSupport, const std::unique_ptr<ScalarFunctionBase_2D> &integrand) const
{
    return integrate(firstSupport, secondSupport, *integrand);
}
