#ifndef INTEGRAL_OPERATOR
#define INTEGRAL_OPERATOR

#include <MyTypes.h>
#include <Operator.h>
#include <mutex>
#include <vector>
#include <memory>


//////////////////////////////////////////
// Base classes for integral operators. //
//////////////////////////////////////////

class DiscreteSpaceOnCurve_1D;
class DiscreteFunction_1D;
class DiscreteFunctionMesh;
class DiscreteSpaceMesh;
class Integrator_1D;
class Integrator_2D;
class Mesh1D;

/**
 * @brief: Base class for integral operator on discrete spaces.
 * DEPRECATED in favor of integralOperatorMesh.
 */
class IntegralOperator : public Operator {
public:
    IntegralOperator(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace);
    virtual ~IntegralOperator();
protected:
    const DiscreteSpaceOnCurve_1D &_trialSpace;
    const DiscreteSpaceOnCurve_1D &_testSpace;
    std::unique_ptr<Integrator_1D> _integrator1D;
    std::unique_ptr<Integrator_2D> _integrator2D;
};

/**
 * @brief: Base class for integral operators on spaces on meshes.
 */
class IntegralOperatorMesh : public Operator {
public:
    IntegralOperatorMesh(const DiscreteSpaceMesh &trialSpace, const DiscreteSpaceMesh &testSpace);
    virtual ~IntegralOperatorMesh();
protected:
    const DiscreteSpaceMesh &_trialSpace;
    const DiscreteSpaceMesh &_testSpace;
    const Mesh1D &_trialMesh;
    const Mesh1D &_testMesh;
    std::unique_ptr<Integrator_1D> _integrator1D;
    std::unique_ptr<Integrator_2D> _integrator2D;
};


#endif
