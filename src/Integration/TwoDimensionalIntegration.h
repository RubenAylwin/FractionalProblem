#ifndef TWO_DIM_INTEGRATION
#define TWO_DIM_INTEGRATION

#include <vector>
#include <memory>
#include <Point.h>
#include <MyTypes.h>
#include <boost/math/quadrature/gauss.hpp>

////////////////////////////////////////////////////////
// Classes for integrating functions of two variables //
////////////////////////////////////////////////////////

//Forward declarations
class ScalarFunctionBase_2D;

/**
 * @brief: Base class for two dimensional integration.
 */
class Integrator_2D {
 public:
    virtual BEM::Complex integrate(const double firstLowerLimit, const double firstUpperLimit, const double secondLowerLimit, const double secondUpperLimit, const ScalarFunctionBase_2D &integrand) const;
    virtual BEM::Complex integrate(const double firstLowerLimit, const double firstUpperLimit, const double secondLowerLimit, const double secondUpperLimit, const std::unique_ptr<ScalarFunctionBase_2D> &integrand) const;
    virtual BEM::Complex integrate(const BEM::Interval1D firstInterval, const BEM::Interval1D secondInterval, const ScalarFunctionBase_2D &integrand) const;
    virtual BEM::Complex integrate(const BEM::Interval1D firstInterval, const BEM::Interval1D secondInterval, const std::unique_ptr<ScalarFunctionBase_2D> &integrand) const;
    virtual BEM::Complex integrate(const BEM::Support1DL firstSupport, const BEM::Support1DL secondSupport, const ScalarFunctionBase_2D &integrand) const;
    virtual BEM::Complex integrate(const BEM::Support1DL firstSupport, const BEM::Support1DL secondSupport, const std::unique_ptr<ScalarFunctionBase_2D> &integrand) const;
    virtual ~Integrator_2D() = default;
 protected:
    Integrator_2D(unsigned n, double firstLowerLimitQuad, double firstUpperLimitQuad, double secondLowerLimitQuad, double secondUpperLimitQuad);
    std::vector<double> _weights;
    std::vector<Point2D> _points;
    double _firstLowerLimitQuad;
    double _firstUpperLimitQuad;
    double _secondLowerLimitQuad;
    double _secondUpperLimitQuad;
};


/**
 * @brief: Gauss-Legendre integration (tensorized). Template class on the number of points for each dim.
 */
template <unsigned N, unsigned M>
class GaussLegendre_2D : public Integrator_2D {
public:
    GaussLegendre_2D(void);
};

/**
 * @brief: Helper class for saving GL points and weights in a more comfortable format.
 */
template <unsigned N, typename C1>
static void emptyGLIntoGivenContainers(C1 &pointContainer, C1 &weightContainer) {
    auto points = boost::math::quadrature::gauss<double, N>::abscissa();
    auto weights = boost::math::quadrature::gauss<double, N>::weights();
    if (N%2 == 0) {
        for (unsigned i = 0; i < N/2; ++i) {
            pointContainer[i] = points[i];
            weightContainer[i] = weights[i];
            pointContainer[i + N/2] = -points[i];
            weightContainer[i + N/2] = weights[i];
        }
    } else {
        pointContainer[0] = points[0];
        weightContainer[0] = weights[0];

        for (unsigned i = 1; i < (N-1)/2 + 1; ++i) {
            pointContainer[i] = points[i];
            weightContainer[i] = weights[i];
            pointContainer[i + (N-1)/2] = -points[i];
            weightContainer[i + (N-1)/2] = weights[i];
        }        
    }

}

/**
 * @brief: Constructor. Only the template parameters are required.
 */
template <unsigned N, unsigned M>
GaussLegendre_2D<N, M>::GaussLegendre_2D(void) :
    Integrator_2D(N*M, -1, 1, -1, 1)
{
    std::vector<double> firstPoints(N, 0), secondPoints(M, 0);
    std::vector<double> firstWeights(N, 0), secondWeights(M, 0);
    emptyGLIntoGivenContainers<N>(firstPoints, firstWeights);
    emptyGLIntoGivenContainers<M>(secondPoints, secondWeights);

    unsigned index = 0;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            _points[index].setX(firstPoints[i]);
            _points[index].setY(secondPoints[j]);
            _weights[index] = firstWeights[i]*secondWeights[j];
            ++index;
        }
    }
}   

#endif
