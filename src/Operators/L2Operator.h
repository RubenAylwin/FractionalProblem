#ifndef L2_OPERATOR
#define L2_OPERATOR

#include <ScalarValuedFunction.h>
#include <DifferentialOperator.h>

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

#endif
