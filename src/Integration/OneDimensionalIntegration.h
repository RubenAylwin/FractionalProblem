#ifndef ONE_DIM_INTEGRATION
#define ONE_DIM_INTEGRATION

///////////////////////////////////////////////////////
// Classes for integrating functions of one variable //
///////////////////////////////////////////////////////

#include <MyTypes.h>
#include <vector>
#pragma push_macro("msg")
#undef msg
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#pragma pop_macro("msg")

// Forward declarations
class ScalarFunctionBase_1D;

/**
 * @brief: Base class for one dimensional integration.
 **/
class Integrator_1D {
 public:
    virtual BEM::Complex integrateSq(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const;
    virtual BEM::Complex integrateSqAnti(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const;
    virtual BEM::Complex integrate(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const;
    virtual BEM::Complex integrate(const double lowerLimit, const double upperLimit, const std::unique_ptr<ScalarFunctionBase_1D> integrand) const;
    virtual BEM::Complex integrate(const BEM::Interval1D support, const ScalarFunctionBase_1D &integrand) const;
    virtual BEM::Complex integrate(const BEM::Interval1D support, const std::unique_ptr<ScalarFunctionBase_1D> integrand) const;
    virtual BEM::Complex integrate(const BEM::Support1DL support, const ScalarFunctionBase_1D &integrand) const;
    virtual BEM::Complex integrate(const BEM::Support1DL support, const std::unique_ptr<ScalarFunctionBase_1D> integrand) const;
    virtual std::vector<double> points(const double lowerLimit, const double upperLimit);
    virtual ~Integrator_1D(void);
 protected:
    Integrator_1D(unsigned n, double lowerLimitQuad, double upperLimitQuad);
    std::vector<double> _weights;
    std::vector<double> _points;
    double _lowerLimitQuad;
    double _upperLimitQuad;
};

/**
 * @brief: Trapezoidal quadrature.
 */
class TrapezoidalQuadrature_1D : public Integrator_1D {
public:
    TrapezoidalQuadrature_1D(unsigned n);
};


/**
 * @brief: Gauss-Legendre quadrature.
 */
class GaussLegendre_1D : public Integrator_1D {
public:
    GaussLegendre_1D(unsigned n);
};


/**
 * @brief: Quadrature for integrands with logaritmic singularities.
 * @TODO: NOT IMPLEMENTED YET.
 */
class LogQuad_1D : public Integrator_1D {
public:
    LogQuad_1D(unsigned n);
};

/**
 * @brief: Adaptive trapezoidal quadrature (wrapper to boost method).
 */
class AdaptiveTrapezoidal_1D : public Integrator_1D {
public:
    AdaptiveTrapezoidal_1D(double tolerance = 1e-06);
    BEM::Complex integrate(const double lowerLimit, const double upperLimit, const ScalarFunctionBase_1D &integrand) const override;
private:
    double _tolerance;
};
#endif
