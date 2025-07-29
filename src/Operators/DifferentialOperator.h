#ifndef DIFFERENTIAL_OPERATOR
#define DIFFERENTIAL_OPERATOR

#include <MyTypes.h>
#include <Operator.h>
#include <mutex>
#include <vector>
#include <memory>

/////////////////////////////////////////////////////
// Class for differential operators (base classes) //
/////////////////////////////////////////////////////

//Forward declarations
class DiscreteSpaceOnCurve_1D;
class DiscreteFunction_1D;
class Integrator_1D;
class Integrator_2D;
class DiscreteSpaceMesh;

/**
 * @brief: Base class. DEPRECATED in favour of DifferentialOperatorMesh.
 */
class DifferentialOperator : public Operator {
public:
    DifferentialOperator(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace);
    virtual ~DifferentialOperator();
protected:
    const DiscreteSpaceOnCurve_1D &_trialSpace;
    const DiscreteSpaceOnCurve_1D &_testSpace;
    std::unique_ptr<Integrator_1D> _integrator1D=nullptr;
    std::unique_ptr<Integrator_2D> _integrator2D=nullptr;
};

/**
 * @brief: Base class.
 */
class DifferentialOperatorMesh : public Operator {
public:
    DifferentialOperatorMesh(const DiscreteSpaceMesh &trialSpace, const DiscreteSpaceMesh &testSpace);
    virtual ~DifferentialOperatorMesh();
protected:
    const DiscreteSpaceMesh &_trialSpace;
    const DiscreteSpaceMesh &_testSpace;
    std::unique_ptr<Integrator_1D> _integrator1D=nullptr;
    std::unique_ptr<Integrator_2D> _integrator2D=nullptr;
};

#endif
