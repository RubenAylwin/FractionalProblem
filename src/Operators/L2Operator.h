#ifndef L2_OPERATOR
#define L2_OPERATOR

#include <MyTypes.h>
#include <ScalarValuedFunction.h>
#include <DifferentialOperator.h>
#include <DiscreteSpaceMesh.h>

class DiscreteSpaceOnCurve_1D;
class ScalarValuedFunction_1D;

class L2 : public DifferentialOperator
{
public:
    L2(ExplicitScalarFunction_1D &function, const DiscreteSpaceOnCurve_1D &trialSpace, DiscreteSpaceOnCurve_1D &testSpace);
    void assembleMassMatrix(void) override;
private:
    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
    ExplicitScalarFunction_1D &_kernel;
};


/**
 * @brief: L2 operator.
 */
class L2Mesh : public DifferentialOperatorMesh {
public:
    enum Side {
        LEFT,
        RIGHT
    };
    L2Mesh(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test);

    template <typename T>
    L2Mesh(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test, T &dFun):
        L2Mesh(trial, test)
    {
        _dFun.reset(new T(dFun));
    };

    L2Mesh(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test, std::shared_ptr<ScalarFunctionBase_1D> &dFun):
        L2Mesh(trial, test)
    {
        _dFun = dFun;
    };

    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
private:
    std::shared_ptr<ScalarFunctionBase_1D> _dFun=nullptr;
};

#endif
