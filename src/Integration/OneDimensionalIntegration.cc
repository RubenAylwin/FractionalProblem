#include <OneDimensionalIntegration.h>
#include <stdexcept>
#include <iostream>
#include <ScalarValuedFunction.h>
#include <cmath>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

/**
 * @brief: Constructor. Number of points, lower and upper limits over which to construct the points.
 */
Integrator_1D::Integrator_1D(unsigned n, double lowerLimitQuad, double upperLimitQuad) :
    _weights(n, 0),
    _points(n, 0),
    _lowerLimitQuad(lowerLimitQuad),
    _upperLimitQuad(upperLimitQuad)
{
}

/**
 * @brief: Destructor.
 */
Integrator_1D::~Integrator_1D(void)
{
}

/**
 * @brief: Integrate the given integrand over the specified limits.
 */
BEM::Complex Integrator_1D::integrate(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const
{
    BEM::Complex result(0, 0);
    const double jacobian = ((upperLimit - lowerLimit)/(_upperLimitQuad - _lowerLimitQuad));;
    for (unsigned i = 0; i < _weights.size(); ++i) {
        double evaluationPoint = lowerLimit + jacobian*(_points[i] - _lowerLimitQuad);
        result += _weights[i]*integrand(evaluationPoint);
    }
    return jacobian*result;
}

/**
 * @brief: Integrate the given integrand over the specified limits with a squared transformation that should help with certain singularities.
 */
BEM::Complex Integrator_1D::integrateSq(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const
{
    auto transformedFunction = [&integrand, lowerLimit](double t) {
        return integrand(t*t - lowerLimit)*2.0*t;
    };
    return integrate(0.0, std::sqrt(upperLimit-lowerLimit), ExplicitScalarFunction_1D(std::move(transformedFunction)));
}

/**
 * @brief: Same as before, but singularity on other endpoint.
 */
BEM::Complex Integrator_1D::integrateSqAnti(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const
{
    auto transformedFunction = [&integrand, upperLimit](double t) {
        return integrand(upperLimit - t*t)*2.0*t;
    };
    return integrate(0.0, std::sqrt(upperLimit-lowerLimit), ExplicitScalarFunction_1D(std::move(transformedFunction)));
}

/**
 * @brief: Integrate over specified interval. Overload for pointer function.
 */
BEM::Complex Integrator_1D::integrate(const double lowerLimit, const double upperLimit, const std::unique_ptr<ScalarFunctionBase_1D> integrand) const
{
    assert(integrand);
    return integrate(lowerLimit, upperLimit, *integrand);
}

/**
 * @brief: Integrate over specified interval. Overload for limits specified with BEM::Interval1D.
 */
BEM::Complex Integrator_1D::integrate(const BEM::Interval1D interval, const ScalarFunctionBase_1D &integrand) const
{
    return integrate(interval.first, interval.second, integrand);
}

/**
 * @brief: Integrate over specified interval. Overload for limits specified with BEM::Interval1D and pointer function.
 */
BEM::Complex Integrator_1D::integrate(const BEM::Interval1D interval, const std::unique_ptr<ScalarFunctionBase_1D> integrand) const
{
    return integrate(interval, *integrand);
}

/**
 * @brief: Overload for BEM::Support1DL (list of intervals over which to integrate).
 */
BEM::Complex Integrator_1D::integrate(const BEM::Support1DL support, const ScalarFunctionBase_1D &integrand) const
{
    BEM::Complex result = 0;
    for (const auto &interval : support) {
        result += integrate(interval, integrand);
    }
    return result;
}

/**
 * @brief: Overload for BEM::Support1DL and pointer function.
 */
BEM::Complex Integrator_1D::integrate(const BEM::Support1DL support, const std::unique_ptr<ScalarFunctionBase_1D> integrand) const
{
    return integrate(support, *integrand);
}

/**
 * @brief: Return the integration points over a specified interval.
 */
std::vector<double> Integrator_1D::points(const double lowerLimit, const double upperLimit)
{
    std::vector<double> result;
    for(size_t i = 0; i < _points.size(); ++i) {
        result.push_back(lowerLimit + (_points[i] - _lowerLimitQuad)*(upperLimit - lowerLimit)/(_upperLimitQuad - _lowerLimitQuad));
    }
    return result;
}


/**
 * @brief: Constructor. Only requires number of points.
 */
TrapezoidalQuadrature_1D::TrapezoidalQuadrature_1D(unsigned n) :
    Integrator_1D(n, 0, 1)
{
    if (n < 2) {
        throw std::invalid_argument("Number of quad. points for trapezoidal rule should be at least 2.");
    }

    for (unsigned i = 0; i < n; ++i) {
        Integrator_1D::_points[i] = i/(n - 1.0);
        Integrator_1D::_weights[i] = 1.0/(n - 1);
    }
    Integrator_1D::_weights[0] *= 0.5;
    Integrator_1D::_weights[n - 1] *= 0.5;
}


/**
 * @brief: Helper templetized class for Gauss-Legendre quadrature.
 * @description: Boost only has the GL quadrature implemented as a template class over integration points.
 * In order to be able to decide integration points on runtime, we have this constructor class that instantiates the given boost GL quadrature and saves the points in a more comfortable format.
 */
template<int N>
class GaussLegendreQuadContainer {
public:
    GaussLegendreQuadContainer() :
        _points(N, 0),
        _weights(N, 0)
    {
        auto points = boost::math::quadrature::gauss<double, N>::abscissa();
        auto weights = boost::math::quadrature::gauss<double, N>::weights();

        if (N%2 == 0) {
            for (unsigned i = 0; i < N/2; ++i) {
                _points[i] = points[i];
                _weights[i] = weights[i];
                _points[i + N/2] = -points[i];
                _weights[i + N/2] = weights[i];
            }
        } else {
            _points[0] = points[0];
            _weights[0] = weights[0];

            for (unsigned i = 1; i < (N-1)/2 + 1; ++i) {
                _points[i] = points[i];
                _weights[i] = weights[i];
                _points[i + (N-1)/2] = -points[i];
                _weights[i + (N-1)/2] = weights[i];
            }        
        }
    }    
    std::vector<double> _points;
    std::vector<double> _weights;    
};

/**
 * @brief: This macro will instantiate the container for various n's, so that then we can just call them later on.
 */
#define QuadInit(z, param, n) {                             \
        if (n == param) {                                   \
            auto Q = GaussLegendreQuadContainer<param>();   \
            for (unsigned i = 0; i < param; ++i) {               \
                _points[i] = Q._points[i];                  \
                _weights[i] = Q._weights[i];                \
            }                                               \
            return;                                         \
        }                                                   \
    }

/**
 * @brief: Constructor. Requires only number of points.
 */
GaussLegendre_1D::GaussLegendre_1D(unsigned n) :
    Integrator_1D(n, -1, 1)
{
    if (n == 0) {
        throw std::invalid_argument("Invalid num of 1D legendre points (zero)");
    }
    BOOST_PP_REPEAT(220, QuadInit, n);
    throw std::invalid_argument("Invalid num of 1D legendre points (max is 200)");
}

/**
 * @brief: Constructor. Requires only number of points.
 * TODO: NOT IMPLEMENTED YET.
 */
LogQuad_1D::LogQuad_1D(unsigned n) :
    Integrator_1D(n, -1, 1)
{
    if (n == 0) {
        throw std::invalid_argument("Invalid num of 1D log points (zero)");
    }
    
    throw std::invalid_argument("Invalid num of 1D log points (not implemented yet)");
}

/**
 * @brief: Constructor. Requires only number of points.
 */
AdaptiveTrapezoidal_1D::AdaptiveTrapezoidal_1D(double tolerance) :
    Integrator_1D(-1, 1, 1),
    _tolerance(tolerance)
{
}

/**
 * @brief: Special integration for adaptive trapezoidal. Wrapper around boost.
 */
BEM::Complex AdaptiveTrapezoidal_1D::integrate(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const {
    return boost::math::quadrature::trapezoidal([&integrand](double t) {return integrand(t);}, lowerLimit, upperLimit, _tolerance, 200);
}
