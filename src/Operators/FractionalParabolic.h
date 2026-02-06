#ifndef FRACTIONAL_PARABOLIC
#define FRACTIONAL_PARABOLIC

#include <MyTypes.h>
#include <DifferentialOperator.h>
#include <vector>
#include <mutex>
#include <memory>

class DiscreteSpaceOnCurve_1D;
class ScalarFunctionBase_1D;
class ExplicitScalarFunction_2D;
class ExplicitScalarFunction_1D;
class Integrator_1D;

/////////////////////////////////////////////////
// ONLY temporal part for fractional problems. //
/////////////////////////////////////////////////

/**
 * @brief: Fractional Parabolic operator
 */
class FractionalParabolic : public DifferentialOperatorMesh {
public:
    FractionalParabolic(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test, int order);

    template <typename T>
    FractionalParabolic(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test, int order, T &dFun):
        FractionalParabolic(trial, test, order)
    {
        _dFun.reset(new T(dFun));
    };

    FractionalParabolic(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test, int order, std::shared_ptr<ScalarFunctionBase_1D> &dFun):
        FractionalParabolic(trial, test, order)
    {
        _dFun = dFun;
    };

    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
private:
    const int _order;
    std::shared_ptr<ScalarFunctionBase_1D> _dFun=nullptr;
};

#endif
