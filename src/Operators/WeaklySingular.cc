#include <WeaklySingular.h>
#include <ScalarValuedFunction.h>
#include <TwoDimensionalIntegration.h>
#include <RegularP0.h>
#include <Mesh.h>
#include <DiscreteSpaceMesh.h>
#include <Green.h>
#include <iostream>
#include <tbb/parallel_for.h>
#include <cassert>
#include <Utilities.h>

/**
 * @brief: Weakly singular (V) operator for the given green function.
 */
WeaklySingular::WeaklySingular(const DiscreteSpaceOnCurve_1D &testSpace, const DiscreteSpaceOnCurve_1D &trialSpace, const GreenFunction2D &greenFunc) :
    IntegralOperator(testSpace, trialSpace),
    _greenFunc{greenFunc}
{
    assert(&(_trialSpace.getCurve()) == &(_testSpace.getCurve()));
    _integrator2D.reset(new GaussLegendre_2D<5, 4>());
}

/**
 * @brief: Explicit destructor
 */
WeaklySingular::~WeaklySingular(void)
{
}

void WeaklySingular::setHighPrecision(void)
{
    _integrator2D.reset(new GaussLegendre_2D<25, 24>());
}

/**
 * @brief: Return const reference to the green function.
 */
const GreenFunction2D &WeaklySingular::getGreen(void)
{
    return _greenFunc;
}

/**
 * Here we wish to compute \int_X \int_Y G(X, Y)*trial(Y)*dY *conj(test(X))*dX
 * We can extract the singular part S(X, Y) so we first numerically compute
 * \int_X \int_Y (G(X, Y) - S(X, Y))*trial(Y)*dY *conj(test(X))*dX however we
 * prefer.
 */
BEM::Complex WeaklySingular::indexedDuality(const unsigned i, const unsigned j)
{
    assert(i < _testSpace.getSize());
    assert(j < _trialSpace.getSize());

    auto &testFunction = _testSpace.basisFunction(i);
    auto &trialFunction = _trialSpace.basisFunction(j);

    auto testSupport = testFunction.brokenSupport();
    auto trialSupport = trialFunction.brokenSupport();
    const auto &curve = _trialSpace.getCurve();
    auto nonSingularIntegrand = ExplicitScalarFunction_2D([&](double t, double s) -> std::complex<double> {
                Point2D auxT(t, 0), auxS(s, 0);
                return (_greenFunc(curve.evaluateAt(t), curve.evaluateAt(s)) - _greenFunc.singularity(auxT, auxS))*curve.jacobian(t)*curve.jacobian(s)*testFunction(t)*trialFunction(s);
        });
    BEM::Complex smoothResult = _integrator2D->integrate(testSupport, trialSupport, nonSingularIntegrand);
    BEM::Complex singularResult = _greenFunc.integrateSingularity(curve, testFunction, trialFunction);
    return smoothResult + singularResult;
}


/**
 * @brief: Weakly singular (V) operator for the given green function.
 */
WeaklySingularMesh::WeaklySingularMesh(const DiscreteSpaceMesh &testSpace, const DiscreteSpaceMesh &trialSpace, const GreenFunction2D &greenFunc) :
    IntegralOperatorMesh(testSpace, trialSpace),
    _greenFunc{greenFunc}
{
    _integrator2D.reset(new GaussLegendre_2D<5, 4>());
}

/**
 * @brief: Weakly singular (V) operator for the given green function.
 */
WeaklySingularMesh::WeaklySingularMesh(const DiscreteSpaceMesh &testSpace, const DiscreteSpaceMesh &trialSpace, const GreenFunction2D &greenFunc, Transformation_2D &&transformation, ExplicitScalarFunction_2D &&jacobian) :
    WeaklySingularMesh(testSpace, trialSpace, greenFunc)
{
    _transformation.reset(new Transformation_2D(transformation));
    _jacobian.reset(new ExplicitScalarFunction_2D(jacobian));                                                                             
}

/**
 * @brief: Explicit destructor
 */
WeaklySingularMesh::~WeaklySingularMesh(void)
{
}

/**
 * @brief: Increase the precision of the integration. Mostly used for testing.
 */
void WeaklySingularMesh::setHighPrecision(void)
{
    _integrator2D.reset(new GaussLegendre_2D<25, 24>());
}

/**
 * @brief: Return const reference to the green function.
 */
const GreenFunction2D &WeaklySingularMesh::getGreen(void)
{
    return _greenFunc;
}

/**
 * Here we wish to compute \int_X \int_Y G(X, Y)*trial(Y)*dY *conj(test(X))*dX
 * We can extract the singular part S(X, Y) so we first numerically compute
 * \int_X \int_Y (G(X, Y) - S(X, Y))*trial(Y)*dY *conj(test(X))*dX however we
 * prefer.
 */
BEM::Complex WeaklySingularMesh::indexedDuality(const unsigned i, const unsigned j)
{
    // TODO: these should throw exceptions.
    assert(i < _testSpace.getSize());
    assert(j < _trialSpace.getSize());

    auto &testFunction = _testSpace.basisFunction(i);
    auto &trialFunction = _trialSpace.basisFunction(j);

    auto testSupportIndices = testFunction.support();
    auto trialSupportIndices = trialFunction.support();

    auto *testElement = &(_testMesh.getElement(0u));
    auto *trialElement = &(_trialMesh.getElement(0u));

    // Right now this only works if given a transformation that parametrizes the curve.
    if (_transformation) {
        auto &curve = *_transformation;
        auto nonSingularIntegrand = ExplicitScalarFunction_2D([&](double t, double s) -> std::complex<double> {
                Point2D pointT = (*testElement)(t);
                Point2D pointS = (*trialElement)(s);
                Point2D auxT((*testElement)[t], 0.0);
                Point2D auxS((*trialElement)[s], 0.0);
                // Right now we remove only the middle singularity, since it appears the others have no real effect.
                return (_greenFunc(curve(pointT), curve(pointS)) - _greenFunc.singularitySimple(auxT, auxS))*testFunction.evaluate(testElement->getIndex(),t)*trialFunction.evaluate(trialElement->getIndex(), s);
            });
        BEM::Complex smoothResult = 0.0;
        BEM::Complex singularResult = 0.0;
        
        // Integrate over the required intervals
        for (unsigned testIndex : testSupportIndices) {
            for (unsigned trialIndex : trialSupportIndices) {
                testElement = &(_testMesh.getElement(testIndex));
                trialElement = &(_trialMesh.getElement(trialIndex));
                smoothResult += _integrator2D->integrate(0, 1, 0, 1, nonSingularIntegrand)*testElement->getSize()*trialElement->getSize();
                singularResult += _greenFunc.integrateSingularity(*testElement, *trialElement, testFunction, trialFunction);
            }
        }
        return smoothResult + singularResult;
    }
    assert(false and "IndexedDuality without transformation not implemented yet.");
    return 0.0;
}

