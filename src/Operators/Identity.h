#ifndef OPERATOR_IDENTITY
#define OPERATOR_IDENTITY
#include <MyTypes.h>
#include <IntegralOperator.h>

////////////////////////////////////////////////////////////
// Identity operator (BEM) for discrete spaces on curves. //
////////////////////////////////////////////////////////////

// Forward declaration.
class DiscreteSpaceOnCurve_1D;

class Identity : public IntegralOperator {
public:
    Identity(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace);
    ~Identity(void);
    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
};

#endif
