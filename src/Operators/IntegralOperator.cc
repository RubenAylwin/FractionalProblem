#include <IntegralOperator.h>
#include <TwoDimensionalIntegration.h>
#include <OneDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <DiscreteSpaceMesh.h>
#include <cassert>

/**
 * @brief: Integral operator on the given spaces (parent class). 
 */
IntegralOperator::IntegralOperator(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace) :
    Operator(testSpace.getSize(), trialSpace.getSize()),
    _trialSpace{trialSpace},
    _testSpace{testSpace},
    _integrator1D{nullptr},
    _integrator2D{nullptr}
{
    assert(&(_trialSpace.getCurve()) == &(_testSpace.getCurve()));
}

/**
 * @brief: Explicit destructor.
 */
IntegralOperator::~IntegralOperator(void)
{
}




/**
 * @brief: Integral operator on the given spaces (parent class). 
 */
IntegralOperatorMesh::IntegralOperatorMesh(const DiscreteSpaceMesh &trialSpace, const DiscreteSpaceMesh &testSpace) :
    Operator(testSpace.getSize(), trialSpace.getSize()),
    _trialSpace{trialSpace},
    _testSpace{testSpace},
    _trialMesh{_trialSpace.getMesh()},
    _testMesh{_testSpace.getMesh()},
    _integrator1D{nullptr},
    _integrator2D{nullptr}
{
}

/**
 * @brief: Explicit destructor.
 */
IntegralOperatorMesh::~IntegralOperatorMesh(void)
{
}
