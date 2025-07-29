#include <GreenQP.h>
#include <GreenH.h>
#include <GreenL.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <cassert>
#include <cmath>
#include <numeric>
#include <iostream>

/**
 * @brief: Window for windowed sum.
 */
static double window(double x, double x0, double x1)
{
    if (std::abs(x) <= x0) {
        return 1.0;
    }
    
    if (std::abs(x) >= x1) {
        return 0.0;
    }
    
    double u = (std::abs(x) - x0)/(x1 - x0);
    return std::exp((2.*std::exp(-1.0/u))/(u - 1.0));
}

/**
 * @brief: Constructor.
 */
GreenQP2D::GreenQP2D(double period) :
    _period(period)
{
    assert(_period > 0 and "Period should be positive.");
}

/**
 * @brief: Explicit destructor.
 */
GreenQP2D::~GreenQP2D(void)
{
}

/**
 * @brief: Singularity depending on the base green function's singularity.
 */
BEM::Complex GreenQP2D::singularity(const Point2D &X, const Point2D &Y) const
{
    return _green->singularity(X, Y)
        + (1./_qp)*_green->singularity(X + Point2D(_period, 0), Y)
        + (_qp)*_green->singularity(X - Point2D(_period, 0), Y);
}

/**
 * @brief: Only central singularity (the most important one).
 */
BEM::Complex GreenQP2D::singularitySimple(const Point2D &X, const Point2D &Y) const
{
    return _green->singularity(X, Y);
}

/**
 * @brief: Sum computing the QP function. Hardcoded 100 terms.
 */
BEM::Complex GreenQP2D::displacedSum(const Point2D &X, const Point2D &Y) const
{
    const int sumTerms = 100;
    BEM::Complex sum(0, 0);
    for (int n = -sumTerms; n <= sumTerms; ++n) {
        sum += (*_green)(X - Point2D(_period*n, 0), Y)*std::pow(_qp, n);
    }
    return sum;
}

/**
 * @brief: Windowed sum.
 */
BEM::Complex GreenQP2D::windowedSum(const Point2D &X, const Point2D &Y) const
{
    Point2D r = X - Y;
    BEM::Complex sum(0, 0);

    double c = 0.35;
    double A = _windowTerms*_period;
    std::vector<BEM::Complex> results(2*_windowTerms + 1, 0);
    for (int n = -_windowTerms; n <= _windowTerms; ++n) {
        const Point2D displacement(_period*n, 0);
        sum += (*_green)(X - Point2D(_period*n, 0), Y)*std::pow(_qp, n)*window(r.getX() - _period*n, c*A, A);
    }
    return sum;
}

/**
 * @brief: Integrating the singularity.
 */
BEM::Complex GreenQP2D::integrateSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction) const
{
    auto middleLog = _green->integrateSingularity(curve, testFunction, trialFunction);

    auto leftLog = _green->integrateShiftedSingularity(curve, testFunction, trialFunction, _period);
    
    auto rightLog = _green->integrateShiftedSingularity(curve, testFunction, trialFunction, -_period);
    
    return middleLog + rightLog*_qp + leftLog/_qp;
}

/**
 * @brief: Same as before, but when given mesh elements. Only middle singularity because it turns out we can ignore the remaining ones.
 */
BEM::Complex GreenQP2D::integrateSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction) const
{
    auto middleLog = _green->integrateSingularity(testElement, trialElement, testFunction, trialFunction);
    return middleLog;
}

/**
 * @brief: Evaluation.
 */
BEM::Complex GreenQP2D::operator()(const Point2D &X, const Point2D &Y) const
{
    return windowedSum(X, Y);
}

/**
 * @brief: Constructor.
 */
GreenHQP2D::GreenHQP2D(double period, double angle, double wavenumber) :
    GreenQP2D(period),
    _angle(angle),
    _wavenumber(wavenumber)
{
    assert(wavenumber > 0 and "Error: wavenumber should be positive");
    _green = std::unique_ptr<PeriodizableGreenFunction2D>(new GreenH2D(wavenumber));
    _greenH = static_cast<GreenH2D*>(_green.get());
    _qp = std::exp(BEM::I*wavenumber*std::sin(angle)*period);
}

/**
 * @brief: Constructor with specific quad points.
 */
GreenHQP2D::GreenHQP2D(double period, double angle, double wavenumber, unsigned quadPoints) :
    GreenHQP2D(period, angle, wavenumber)
{
    _greenH->set1DQuad(quadPoints);
}

/**
 * @brief: Increase the precision of integration.
 */
void GreenHQP2D::setHighPrecision(int i)
{
    _greenH->setHighPrecision(i);
}

/**
 * @brief: Other representation of the green function.
 */
BEM::Complex GreenHQP2D::spectralSum(const Point2D &X, const Point2D &Y) const
{

    Point2D r = X - Y;

    const int sumTerms = _windowTerms;
 
    const double K = _wavenumber;
    
    BEM::Complex sum(0, 0);
    for (int n = -sumTerms; n <= sumTerms; ++n) {
        double Bn = K*std::sin(_angle) + n*2.0*M_PI/(_period);
        BEM::Complex Gn = (K*K >= Bn*Bn) ? std::sqrt(K*K - Bn*Bn) : BEM::I*std::sqrt(Bn*Bn - K*K);
        BEM::Complex xExp = std::exp(BEM::I*Bn*(r.getX()));
        BEM::Complex yExp = std::exp(BEM::I*Gn*std::abs<double>(r.getY()));
        sum += xExp*yExp/Gn;
    }
    return -(BEM::I/(2.*_period))*sum;
}

/**
 * @brief: Constructor.
 */
GreenLQP2D::GreenLQP2D(double period) :
    GreenQP2D(period)
{
    _green = std::unique_ptr<PeriodizableGreenFunction2D>(new GreenL2D());
}
