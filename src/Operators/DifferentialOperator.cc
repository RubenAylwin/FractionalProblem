#include <DifferentialOperator.h>
#include <TwoDimensionalIntegration.h>
#include <OneDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <DiscreteSpaceMesh.h>
#include <cassert>
#include <memory>

/**
 * @brief: Differential operator on the given spaces (parent class). 
 */
DifferentialOperator::DifferentialOperator(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace) :
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
DifferentialOperator::~DifferentialOperator(void)
{
}


/**
 * @brief: Differential operator on the given spaces (parent class). 
 */
DifferentialOperatorMesh::DifferentialOperatorMesh(const DiscreteSpaceMesh &trialSpace, const DiscreteSpaceMesh &testSpace) :
    Operator(testSpace.getSize(), trialSpace.getSize()),
    _trialSpace{trialSpace},
    _testSpace{testSpace},
    _integrator1D{nullptr},
    _integrator2D{nullptr}
{
}

/**
 * @brief: Explicit destructor.
 */
DifferentialOperatorMesh::~DifferentialOperatorMesh(void)
{
}
