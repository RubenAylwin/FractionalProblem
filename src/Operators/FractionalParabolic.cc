#include <FractionalParabolic.h>
#include <TwoDimensionalIntegration.h>
#include <OneDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <ScalarValuedFunction.h>
#include <iostream>
#include <ranges>
#include <Utilities.h>
#include <Msg.h>
#include <set>

useMessages("FRAC_PAR");
/**
 * @brief: Constructor.
 */
FractionalParabolic::FractionalParabolic(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test, int order):
    DifferentialOperatorMesh(trial, test),
    _order{order}
{
    _integrator1D.reset(new GaussLegendre_1D(5));//Hardcoded. Should be enough (maybe even 2 in enough).
}

BEM::Complex FractionalParabolic::indexedDuality(const unsigned i, const unsigned j)
{
    auto &test = _testSpace.basisFunction(i);
    auto &trial = _trialSpace.basisFunction(j);
    auto &trialL = trial.derivative(_order);
    auto &mesh = _testSpace.getMesh();
    BEM::Complex result = 0.0;
    auto integrationElementSize = mesh.getElement(0).getSize();
    unsigned l = test.support()[0];
    ExplicitScalarFunction_1D integrand([&](double t) {
        auto point = mesh.getElement(l)(t);
        return test.evaluate(l, t)*trialL.evaluate(l, t)*integrationElementSize*(*_dFun)(point.getX());
    }
        );

    // Right now we integrate over all elements, but we know that the integration domain is less that this.
    for (unsigned k : test.support()) {
        l=k;
        integrationElementSize = mesh.getElement(l).getSize();
        result += _integrator1D->integrate(0, 1, integrand);
    }
    
    return result;
}

