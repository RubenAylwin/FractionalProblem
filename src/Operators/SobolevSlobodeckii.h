#ifndef SOBOLEV_SLOBODECKII
#define SOBOLEV_SLOBODECKII

#include <MyTypes.h>
#include <DifferentialOperator.h>

///////////////////////////////////////////
// Class for Sobolev Slobodeckii product //
///////////////////////////////////////////


// Forward declaration
class DiscreteSpaceOnCurve_1D;

class SobolevSlobodeckii : public DifferentialOperator {
public:
    SobolevSlobodeckii(const DiscreteSpaceOnCurve_1D &space, double order);
    BEM::Complex indexedDuality(const unsigned i, const unsigned j);
private:
    const double _order;
    BEM::Support1DL _elements;
};

#endif
