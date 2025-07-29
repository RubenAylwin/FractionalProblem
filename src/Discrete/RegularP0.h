#ifndef REGULAR_P0
#define REGULAR_P0

#include <Curve.h>
#include <utility>
#include <ScalarValuedFunction.h>
#include <Eigen/Dense>
#include <vector>
#include <MyTypes.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <mutex>
#include <cassert>
#include <DiscreteSpace.h>
// Forward declarations
class DiscreteSpace_1D;
class Integrator_1D;
class Integrator_2D;


//////////////////////////////////////////////////////////////////
// Classes for spaces of piecewise P0 functions on a partition. //
/////////////////////////////////////////////////////////////////


/**
 * @brief: Space of piecewise P0 functions on a curve and partition.
 * DEPRECATED
 */
class RegularP0_1D : public DiscreteSpaceOnCurve_1D {
public:
    RegularP0_1D(PolyPeriodicCurve &curve, int size);
    RegularP0_1D(const Curve2D &curve, int size);
    ~RegularP0_1D(void);
    unsigned getSize(void) const override { return _size; };
    DiscreteFunctionPtr generateFunction(const std::vector<BEM::Complex> coefficients) const override;
    const BasisFunction_1D &basisFunction(const int globalNumber) const override;
    const DiscreteFunctional_1D &degreeOfFreedom(const int globalNumber) const override;
    const std::string &getBaseName(void) const override { return _baseName; };
    const std::vector<double> &getPartition(void) const override {return _partition;};
    unsigned globalIndex(unsigned element, unsigned localFunction) const override;
    // Basis
    class P0Function : virtual public DiscreteFunction_1D {
    public:
        BEM::Complex operator()(const double t) const override;
        BEM::Interval1D support(void) const override;
        const BEM::Support1DL &brokenSupport(void) const override;
    private:
        P0Function(const std::vector<double> &partition, const std::vector<BEM::Complex> coefficients, const Curve2D &curve);
        const std::vector<double> &_partition;
        const std::vector<BEM::Complex> _coefficients;
        const Curve2D &_curve;
        friend class RegularP0_1D;
        mutable std::mutex _derMutex;
        mutable std::mutex _supportMutex;
    };

    class BasisFunction : public BasisFunction_1D, public P0Function {
    private:
        BasisFunction(const std::vector<double> &partition, int index, const Curve2D &curve);
        ScalarFunctionBase_2D* antiDerivative(Direction direction) const;
        ScalarFunctionBase_1D *derivative(void) const;
        mutable std::unique_ptr<ScalarFunctionBase_2D> _antiDerivativeLeft = nullptr;
        mutable std::unique_ptr<ScalarFunctionBase_2D> _antiDerivativeRight = nullptr;
        mutable std::unique_ptr<ScalarFunctionBase_1D> _derivative = nullptr;

        friend class RegularP0_1D;
    };

private:
    const int _size;
    std::vector<double> _partition;
    const std::string _baseName = "RegularP0";
    mutable std::unordered_map<int, std::unique_ptr<BasisFunction_1D>> _basis;
    mutable std::unordered_map<int, std::unique_ptr<DiscreteFunctional_1D>> _dof;
    mutable std::mutex _basisMutex;
    mutable std::mutex _dofMutex;
};

#endif
