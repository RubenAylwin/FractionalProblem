#include <GreenLog.h>
#include <Point.h>
#include <DiscreteSpace.h>
#include <DiscreteSpaceMesh.h>
#include <Curve.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <Mesh.h>
#include <iostream>
#include <Msg.h>
useMessages("GREEN_LOG");

/**
 * @brief: Constructor.
 */
GreenLogSing2D::GreenLogSing2D(void) :
    _integrator1D(new GaussLegendre_1D(8)),
    _integrator2D(new GaussLegendre_2D<5,5>())
{
}

/**
 * @brief: Explicit destructor.
 */
GreenLogSing2D::~GreenLogSing2D(void)
{
}

/**
 * @brief: Set points of one dimensional quadrature.
 */
void GreenLogSing2D::set1DQuad(unsigned N)
{
    _integrator1D.reset(new GaussLegendre_1D(N));
}

/**
 * @brief: Set quadratures to high(er) precision.
 */
void GreenLogSing2D::setHighPrecision(int level)
{
    if (level == 0) {
        _integrator1D.reset(new GaussLegendre_1D(40));
        _integrator2D.reset(new GaussLegendre_2D<25, 24>());
    } if (level == 1) {
        _integrator1D.reset(new GaussLegendre_1D(80));
        _integrator2D.reset(new GaussLegendre_2D<50, 51>());
    } if (level == 2) {
        _integrator1D.reset(new GaussLegendre_1D(200));
        _integrator2D.reset(new GaussLegendre_2D<401, 400>());
    }
}

/**
 * @brief: Singularity to extract.
 */
BEM::Complex GreenLogSing2D::singularity(const Point2D &X, const Point2D &Y) const
{
    return (1.0/(2.0*M_PI))*std::log(std::abs(X.getX() - Y.getX()));
}

/**
 * @brief: Singularity to extract.
 */
BEM::Complex GreenLogSing2D::singularitySimple(const Point2D &X, const Point2D &Y) const
{
    return (1.0/(2.0*M_PI))*std::log(std::abs(X.getX() - Y.getX()));
}

/**
 * @brief: Integration of the Singularity(X,Y)*Trial(Y)*Test(X)
 */
BEM::Complex GreenLogSing2D::integrateSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction) const
{
    msg(5) << "start::GreenLogSing2D::integrateSingularity" << endMsg;
    return integrateShiftedSingularity(curve, testFunction, trialFunction, 0);
    msg(5) << "end::GreenLogSing2D::integrateSingularity" << endMsg;
}

/**
 * @brief: Integration of the Singularity(X,Y)*Trial(Y)*Test(X)
 */
BEM::Complex GreenLogSing2D::integrateSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction) const
{
    msg(5) << "start::GreenLogSing2D::integrateSingularity" << endMsg;
    return integrateShiftedSingularity(testElement, trialElement, testFunction, trialFunction, 0);
    msg(5) << "end::GreenLogSing2D::integrateSingularity" << endMsg;
}

/**
 * @brief: For use by the QP Green Function. Integration of the singularity shifted in the first coordinate.
 */
BEM::Complex GreenLogSing2D::integrateShiftedSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction, const double shift) const
{
    msg(5) << "start::GreenLogSing2D::integrateShiftedSingularity" << endMsg;
    // We know that the singularity can be extracted as a function of the integration on the parametrized space, namely log(|t - s|), so the given curve is not
    // necessary here and we ignore it. We just need to integrate log(|t - s|)*trial(s)*test(t).

    auto trialAntiDerivativeLPtr = trialFunction.antiDerivative(BasisFunction_1D::Direction::LEFT);
    auto trialAntiDerivativeRPtr = trialFunction.antiDerivative(BasisFunction_1D::Direction::RIGHT);
    msg(5) << "Got the antiDerivatives" << endMsg;
    auto &trialAntiDerivativeL = *trialAntiDerivativeLPtr;
    auto &trialAntiDerivativeR = *trialAntiDerivativeRPtr;
    msg(5) << "Got the antiDerivative pointers" << endMsg;
    auto trialSupport = trialFunction.brokenSupport();
    msg(5) << "Got the trial support" << endMsg;
    auto shiftedTestSupport = testFunction.brokenSupport() + shift;
    msg(5) << "Got the shifted test support" << endMsg;

    BEM::Complex oneDimRes = 0.0;
    for (const auto trialSupportInterval : trialSupport) {
        ExplicitScalarFunction_1D oneDimIntegrand([&](double t) -> BEM::Complex {
            if (t == trialSupportInterval.first or t == trialSupportInterval.second) {
                return 0;
            }
            return (std::log(std::abs(t - trialSupportInterval.second))*trialAntiDerivativeL(trialSupportInterval.second, -(t))
                    - std::log(std::abs(t - trialSupportInterval.first))*trialAntiDerivativeR(trialSupportInterval.first, -(t))
                    + std::log(std::abs(t - 0.5*(trialSupportInterval.first+trialSupportInterval.second)))*(trialAntiDerivativeL(0.5*(trialSupportInterval.first+trialSupportInterval.second)-0.00000001, -(t)) - trialAntiDerivativeR(0.5*(trialSupportInterval.first+trialSupportInterval.second)+0.00000001, -(t)))
                    )*testFunction(t - shift)*curve.jacobian(t - shift);
        });
        oneDimRes += _integrator1D->integrate(shiftedTestSupport, oneDimIntegrand);
    }
    
    ExplicitScalarFunction_2D twoDimIntegrand([&](double t, double s) -> BEM::Complex {
        if (s == t) {
            return 0.0;
        }
        auto val = trialAntiDerivativeL(s, -t)*testFunction(t - shift)/(t - s)*curve.jacobian(t - shift);
        return val;
        
    });
    BEM::Complex secondDimRes =_integrator2D->integrate(shiftedTestSupport, trialSupport, twoDimIntegrand);

    msg(5) << "OneDimRes: " << oneDimRes << endMsg;
    msg(5) << "SecDimRes: " << secondDimRes << endMsg;
    msg(5) << "end::GreenLogSing2D::integrateShiftedSingularity" << endMsg;
    
    return (oneDimRes + secondDimRes)/(2.*M_PI);
}

/**
 * @brief: For use by the QP Green Function. Integration of the singularity shifted in the first coordinate.
 */
BEM::Complex GreenLogSing2D::integrateShiftedSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction, const double shift) const
{
    msg(5) << "start::GreenLogSing2D::integrateShiftedSingularity" << endMsg;
    // We know that the singularity can be extracted as a function of the integration on the parametrized space, namely log(|t - s|),
    // so the given curve is not necessary here and we ignore it. We just need to integrate log(|t - s|)*trial(s)*test(t).
    Point2D shiftedPoint(shift, 0.);
    double singleShift = shift == 0.0 ? 0.0 : shift/std::abs(shift);
    assert(singleShift == 0.0);
    const auto &trialAntiDerivative = trialFunction.antiDerivative(trialElement.getIndex());
    double c = trialElement[0];
    double d = trialElement[1];

    ExplicitScalarFunction_1D oneDimIntegrand([&](double t) -> BEM::Complex {
        return (std::log(std::abs(testElement[t] - d))*(trialAntiDerivative(1) - t) - std::log(std::abs(testElement[t] - c))*(trialAntiDerivative(0) - t))*testFunction.evaluate(testElement.getIndex(), -t);
    });
    ExplicitScalarFunction_2D twoDimRem([&](double t, double s) -> BEM::Complex {
        if (s == t) {
            return 0.0;
        }
        return (d-c)*(trialAntiDerivative(s) - t)*testFunction.evaluate(testElement.getIndex(), -t)/(testElement[t] - trialElement[s]);
    });    
    BEM::Complex oneDim = _integrator1D->integrate(0, 1, oneDimIntegrand);
    BEM::Complex secondRem = _integrator2D->integrate(0, 1, 0, 1, twoDimRem);
    return testElement.getSize()*trialElement.getSize()*(oneDim + secondRem)/(2.*M_PI);
}
