#include <FractionalDerivative.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <ScalarValuedFunction.h>
#include <Msg.h>
#include <Utilities.h>

useMessages("FRAC_DER");

/**
 * @brief: Constructor.
 */
LeftFracDerivative::LeftFracDerivative(const DiscreteSpaceMesh &space, int order):
    DifferentialOperatorMesh(space, space),
    _order{order}
{
    _integrator1D.reset(new GaussLegendre_1D(2));
}

/**
 * @brief: L2-product of left derivatives of test function i and trial function j.
 * @desc: This should go in position (i,j) of the associated matrix.
 */
BEM::Complex LeftFracDerivative::indexedDuality(const unsigned i, const unsigned j)
{
    auto &testL = _testSpace.basisFunction(i).derivative(_order);
    auto &trialL = _trialSpace.basisFunction(j).derivative(_order);
    auto &mesh = _testSpace.getMesh();
    BEM::Complex result = 0.0;
    auto integrationElementSize = mesh.getElement(0).getSize();//Just for initialization.
    unsigned l = 0;
    ExplicitScalarFunction_1D integrand([&](double t) {
        return testL.evaluate(l,t)*trialL.evaluate(l,t)*integrationElementSize;
    });

    //TODO: Use the support of the functions.
    for (l = 0 ; l < mesh.numElements(); ++l) {
        integrationElementSize = mesh.getElement(l).getSize();
        result += _integrator1D->integrate(0, 1, integrand);
    }
    return result;
}
