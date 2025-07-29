#include <Identity.h>
#include <ScalarValuedFunction.h>
#include <OneDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <Utilities.h>
#include <cassert>
#include <tbb/parallel_for.h>

/**
 * @brief: Identity operator.
 */
Identity::Identity(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace) :
    IntegralOperator(trialSpace, testSpace)
{
    assert(&(_trialSpace.getCurve()) == &(_testSpace.getCurve()));
    _integrator1D.reset(new GaussLegendre_1D(3));
    
}

/**
 * @brief: Explicit destructor.
 */
Identity::~Identity(void)
{
}

/**
 * @brief: Integration of test(i)*trial(j).
 */
BEM::Complex Identity::indexedDuality(const unsigned i, const unsigned j)
{
    assert(static_cast<size_t>(i) < _testSpace.getSize());
    assert(static_cast<size_t>(j) < _trialSpace.getSize());
    
    auto &trialFunction = _trialSpace.basisFunction(j);
    auto &testFunction = _testSpace.basisFunction(i);
    auto supTrial = trialFunction.support();
    auto supTest = testFunction.support();
    auto sup = BEM::intersect(supTrial, supTest);
    if (sup.first > sup.second) {
        return 0.0;
    }
    auto integrand = (trialFunction)*(testFunction);
    return _integrator1D->integrate(testFunction.brokenSupport(), (*integrand));
}

