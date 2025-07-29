#ifndef DISCRETE_SPACE
#define DISCRETE_SPACE

#include <Curve.h>
#include <Element.h>
#include <ScalarValuedFunction.h>
#include <MyTypes.h>
#include <utility>
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <mutex>
#include <cassert>


///////////////////////////////////////////////////////////////////////
// Classes for discrete spaces and functions.                        //
// Right now this code is DEPRECATED in favour of DiscreteSpaceMesh. //
// This spaces assume a disc. over [0,1] which was enough for the    //
// first application that this code was used for.                    //
///////////////////////////////////////////////////////////////////////


// Forward declarations
class DiscreteSpace_1D;
class Integrator_1D;
class Integrator_2D;

/**
 * @brief: Base class for discrete functions of one variable.
 */
class DiscreteFunction_1D : public ScalarFunctionBase_1D {
public:
    DiscreteFunction_1D(const DiscreteFunction_1D &other) = delete;
    DiscreteFunction_1D(void) {};
    virtual ~DiscreteFunction_1D(void) = default;
    protected:
    mutable std::unique_ptr<BEM::Interval1D> _support;
    mutable std::unique_ptr<BEM::Support1DL> _brokenSupport;
};

/**
 * @brief: Base class for basis functions of discrete spaces.
 */
class BasisFunction_1D : virtual public DiscreteFunction_1D
{
public:
    //Enum
    enum Direction {
        LEFT,
        RIGHT
    };
    BasisFunction_1D(void) = default;
    virtual ScalarFunctionBase_2D *antiDerivative(Direction direction) const = 0;
    virtual ScalarFunctionBase_1D *derivative(void) const = 0;
};

//Typenames
using DiscreteFunctionPtr = std::unique_ptr<DiscreteFunction_1D>;
using DiscreteFunctional_1D = std::function<BEM::Complex (const ScalarFunctionBase_1D &function)>;

/**
 * @brief: Base class for discrete space of one dimensional functions.
 */
class DiscreteSpace_1D {
public:
    DiscreteSpace_1D(void);
    virtual ~DiscreteSpace_1D(void);
    virtual BEM::ColVector testAgainstBasis(const ScalarFunctionBase_1D &function) const;
    virtual BEM::Complex testAgainstBaseElement(const ScalarFunctionBase_1D &function, int globalNumber) const;
    virtual BEM::Matrix identityOp(void) const;
    virtual DiscreteFunctionPtr generateFunction(const BEM::ColVector &coefficients) const;
    virtual BEM::ColVector projectFunction(const ScalarFunctionBase_1D &function) const;
    virtual DiscreteFunctionPtr generateFunction(const std::vector<BEM::Complex> coefficients) const = 0;
    virtual const BasisFunction_1D &basisFunction(const int globalNumber) const = 0;
    virtual const DiscreteFunctional_1D &degreeOfFreedom(const int globalNumber) const = 0;
    virtual const std::vector<double> &getPartition(void) const = 0;
    virtual unsigned getSize(void) const = 0;
    virtual const std::string &getBaseName(void) const = 0;
    virtual unsigned globalIndex(unsigned element, unsigned localFunction) const {return 0*(element + localFunction);};
    virtual const std::vector<Element_1D> &getElements(void) const {return _elements;};
protected:
    std::unique_ptr<Integrator_1D> _integrator1D;
    std::unique_ptr<Integrator_2D> _integrator2D;
    std::vector<Element_1D> _elements;
};

/**
 * @brief: Discrete space over a given curve in 2D space.
 */
class DiscreteSpaceOnCurve_1D : public DiscreteSpace_1D {
public:
    DiscreteSpaceOnCurve_1D(const Curve2D &curve) : DiscreteSpace_1D(), _curve{curve} {}
    const Curve2D &getCurve(void) const { return _curve; }
protected:
    const Curve2D &_curve;
};

/**
 * @brief: P1 space over a curve.
 * DEPRECATED.
 */
class RegularP1_1D : public DiscreteSpaceOnCurve_1D {
public:
    RegularP1_1D(PolyPeriodicCurve &curve, int size, bool zeroBC = false);
    RegularP1_1D(const Curve2D &curve, int size, bool zeroBC = false);
    ~RegularP1_1D(void);
    void setZeroBC(void);
    unsigned getSize(void) const override { return _size; };
    DiscreteFunctionPtr generateFunction(const std::vector<BEM::Complex> coefficients) const override;
    const BasisFunction_1D &basisFunction(const int globalNumber) const override;
    const DiscreteFunctional_1D &degreeOfFreedom(const int globalNumber) const override;
    const std::string &getBaseName(void) const override { return _baseName; };
    const std::vector<double> &getPartition(void) const override {return _partition;};
    unsigned globalIndex(unsigned element, unsigned localFunction) const override;
    
    // Basis
    class P1Function : virtual public DiscreteFunction_1D {
    public:
        BEM::Complex operator()(const double t) const override;
        BEM::Interval1D support(void) const override;
        const BEM::Support1DL &brokenSupport(void) const override;
    private:
        P1Function(const std::vector<double> &partition, const std::vector<BEM::Complex> coefficients, const Curve2D &curve);

        const std::vector<double> &_partition;
        std::vector<BEM::Complex> _coefficients;
        const Curve2D &_curve;
        friend class RegularP1_1D;
        mutable std::mutex _derMutex;
    };

    class BasisFunction : public BasisFunction_1D, public P1Function {
    private:
        BasisFunction(const std::vector<double> &partition, unsigned index, const Curve2D &curve);
        ScalarFunctionBase_2D *antiDerivative(Direction direction) const override;
        ScalarFunctionBase_1D *derivative(void) const override;
        unsigned _index = 0;
        mutable std::unique_ptr<ScalarFunctionBase_2D> _antiDerivativeLeft = nullptr;
        mutable std::unique_ptr<ScalarFunctionBase_2D> _antiDerivativeRight = nullptr;
        mutable std::unique_ptr<ScalarFunctionBase_1D> _derivative = nullptr;

        friend class RegularP1_1D;
    };

protected:
    const unsigned _size;
    const bool _zeroBC;
    std::vector<double> _partition;
    const std::string _baseName = "RegularP1";
    mutable std::unordered_map<int, std::unique_ptr<BasisFunction_1D>> _basis;
    mutable std::unordered_map<int, std::unique_ptr<DiscreteFunctional_1D>> _dof;
    mutable std::mutex _basisMutex;
    mutable std::mutex _dofMutex;
};

/**
 * @brief: Graded P1 space over a curve.
 * DEPRECATED.
 */
class GeoP1_1D : public RegularP1_1D
{
public:
    GeoP1_1D(PolyPeriodicCurve &curve, int size, double power, bool zeroBC = false);
    GeoP1_1D(const Curve2D &curve, int size, double power, bool zeroBC = false);
private:
    double _power;
};

#endif
