#ifndef FRACTIONAL_DERIVATIVE
#define FRACTIONAL_DERIVATIVE

#include <MyTypes.h>
#include <DifferentialOperator.h>
#include <vector>
#include <mutex>
#include <memory>

///////////////////////////////////////////////////////
// Product of left RL derivatives of the same order. //
///////////////////////////////////////////////////////

class DiscreteSpaceMesh;

class LeftFracDerivative : public DifferentialOperatorMesh {
public:
    LeftFracDerivative(const DiscreteSpaceMesh &space, int order);
    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
private:
    const int _order;
};

#endif
