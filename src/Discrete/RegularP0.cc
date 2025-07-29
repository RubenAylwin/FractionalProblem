#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <ScalarValuedFunction.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <MyTypes.h>
#include <Utilities.h>
#include <Msg.h>
#include <RegularP0.h>

useMessages("REG_P0");


/**
 * @brief: Constructor on a periodic curve.
 */
RegularP0_1D::RegularP0_1D(PolyPeriodicCurve &curve, int size) :
    RegularP0_1D(static_cast<Curve2D &>(curve), size)
{
    auto intrinsicPartition = curve.getPartition();
    assert(intrinsicPartition.size() < _size + 1u);

    for (const auto &val : intrinsicPartition) {
        if (val == 0 or val == 1) {
            continue;
        }
        double integerPart;
        double fractionalPart = std::modf(val*_size, &integerPart);

        if (fractionalPart == 0) {
            continue;
        }
        int index = (fractionalPart > 0.5) ? integerPart + 1 : integerPart;
        _partition[index] = val;
    }
    
    // Correctness checks
    assert(_partition[0] == 0.0);
    assert(_partition[_size] == 1.0);
}


/**
 * @brief: Constructor on a general curve
 */
RegularP0_1D::RegularP0_1D(const Curve2D &curve, int size) :
    DiscreteSpaceOnCurve_1D(curve),
    _size(size),
    _partition(_size + 1, 0.0)
{
    _integrator1D.reset(new GaussLegendre_1D(8));
    _integrator2D.reset(new GaussLegendre_2D<5, 5>());
    // Size of each interval
    double intervalSize = 1.0/_size;
    _partition[0] = 0.;
    
    // Fill in the nodes of the partition
    for (int i = 1; i < _size + 1; ++i) {
        _partition[i] = i*intervalSize;
        _elements.emplace_back(i-1, i, i-1, 1);
    }

    // Correctness checks
    assert(_partition[0] == 0.0);
    assert(_partition[_size] == 1.0);
}

/**
 * @brief: Destructor
 */
RegularP0_1D::~RegularP0_1D(void)
{
}

/**
 * @brief: Returns the i-th basis function.
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const BasisFunction_1D &RegularP0_1D::basisFunction(const int globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber < _size and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid basis index (out of bounds)");
    }

    std::lock_guard<std::mutex> lock(_basisMutex);
    if (_basis.find(globalNumber) == _basis.end()) {
        // Fill in coefficient vector
        std::vector<BEM::Complex> coefficients(_size, std::complex(0.0, 0.0));

        // Coef = 1 for the corresponding basis
        coefficients[globalNumber] = std::complex(1.0, 0.0);

        // Return unique ptr
        _basis.insert({globalNumber ,std::unique_ptr<BasisFunction_1D>(new BasisFunction(_partition, globalNumber, _curve))});
    }
    return *(_basis.at(globalNumber));
}

/**
 * @brief: Returns the i-th degree of freedom..
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const DiscreteFunctional_1D &RegularP0_1D::degreeOfFreedom(const int globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber < _size and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid DoF index (out of bounds)");
    }

    std::lock_guard<std::mutex> lock(_dofMutex); //Expected to be called in parallel, so lock it.
    if (_dof.find(globalNumber) == _dof.end()) {
        double size = _partition.at(globalNumber+1) - _partition.at(globalNumber);
 _dof.insert({globalNumber ,std::unique_ptr<DiscreteFunctional_1D>(new DiscreteFunctional_1D([this, size, globalNumber](const ScalarFunctionBase_1D &function) -> BEM::Complex {
     auto support = BEM::intersect(function.support(), BEM::Interval1D(_partition.at(globalNumber), _partition.at(globalNumber + 1)));
     return (_integrator1D->integrate(support, function*(*_curve.surfaceMeasure())))/size;
 }))});
    }
    return *(_dof.at(globalNumber));
}


/**
 * @brief: Given a coefficient vector return the associated discrete function.
 * This is the only way to generate a discretization/construct a discrete function.
 */
std::unique_ptr<DiscreteFunction_1D> RegularP0_1D::generateFunction(const std::vector<BEM::Complex> coefficients) const
{
    return std::unique_ptr<DiscreteFunction_1D>(new P0Function(_partition, coefficients, _curve));
}

/**
 * @brief: Construct a discrete function from given coefficients.
 */
RegularP0_1D::P0Function::P0Function(const std::vector<double> &partition, const std::vector<BEM::Complex> coefficients, const Curve2D &curve) :
    _partition(partition),
    _coefficients(coefficients),
    _curve(curve)
{
    assert(_coefficients.size() == _partition.size() - 1);
}

/**
 * @brief: Evaluate the function.
 * @Assumption: Parametrization of the curve on [0,1]
 */
BEM::Complex RegularP0_1D::P0Function::operator()(const double t) const
{
    // Check for correct range (TODO: THIS SHOULD BE A THROW)
    assert(t >= 0 and t <= 1);

    // Find the index
    auto index = BEM::findIndex(_partition, t, BasisFunction_1D::Direction::LEFT);
    
    // Return
    return _coefficients[index]/_curve.jacobian(t);
}

/**
 * @brief: Return the support of a given function as an interval.
 */
BEM::Interval1D RegularP0_1D::P0Function::support(void) const
{
    if (_support) {
        return *_support;
    }
    
    // Find first non zero coefficient
    unsigned i = 0;
    while (i <_coefficients.size() and std::abs(_coefficients[i]) < 1E-10) {
        ++i;
    }

    // If got to the end of the list, support is null.
    if (i == _coefficients.size()) {
        return {0,0};
    }

    // Correctness check
    assert(i < _coefficients.size());
    unsigned lower = i;

    // Find the last non zero coefficient
    i = _coefficients.size() - 1;
    while (i > lower and std::abs(_coefficients[i]) < 1E-10) {
        --i;
    }
    int upper = i + 1;

    // return values
    _support.reset(new BEM::Interval1D(_partition.at(lower), _partition.at(upper)));
    return *_support;
}

/**
 * @brief: given an element and a local numbering of the element, get the global numbering of the function.
 */
unsigned RegularP0_1D::globalIndex(unsigned element, unsigned localFunction) const
{
    assert(localFunction == 0); // Only one function per element
    return element;
}

/**
 * @brief: Return the support of a given function as a list of intervals
 */
const BEM::Support1DL &RegularP0_1D::P0Function::brokenSupport(void) const
{
    std::lock_guard<std::mutex> lock(_supportMutex);
    if (_brokenSupport) {
        return *_brokenSupport;
    }
    _brokenSupport.reset(new BEM::Support1DL{support()});
    return *_brokenSupport;
}

/**
 * @brief: return a function evaluating the derivative of the corresponding function. 
 */
ScalarFunctionBase_1D *RegularP0_1D::BasisFunction::derivative(void) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (not _derivative) {
        auto derivative = [](double t) -> BEM::Complex {
            return t*0.0;
        };
        _derivative.reset(new ExplicitScalarFunction_1D(derivative));
    }
    return _derivative.get();
}

/**
 * @brief: return a function evaluating the anti-derivative of the corresponding function. First coordinate corresponds
 * to the point of evaluation. The second coordinate gives the constant (since derivative of constant is 0).
 */
ScalarFunctionBase_2D *RegularP0_1D::BasisFunction::antiDerivative(Direction direction) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (direction == Direction::LEFT) {
        if (not _antiDerivativeLeft) {
            auto antiDerivativeLeft = [&] (double t, double s) -> BEM::Complex {
                // Check for correct range (TODO: THIS SHOULD BE A THROW) rate
                assert(t >= 0 and t <= 1);
                
                auto index = BEM::findIndex(_partition, t, Direction::LEFT);
                return _coefficients[index]*(t + s);
            };
            _antiDerivativeLeft.reset(new ExplicitScalarFunction_2D(antiDerivativeLeft));
        }
        return _antiDerivativeLeft.get();
    } else /*if (direction == Direction::RIGHT)*/ {
        if (not _antiDerivativeRight) {
            auto antiDerivativeRight = [&] (double t, double s) -> BEM::Complex {
                // Check for correct range (TODO: THIS SHOULD BE A THROW) rate
                assert(t >= 0 and t <= 1);
                
                auto index = BEM::findIndex(_partition, t, Direction::RIGHT);
                return _coefficients[index]*(t + s);
            };
            _antiDerivativeRight.reset(new ExplicitScalarFunction_2D(antiDerivativeRight));
        }
        return _antiDerivativeRight.get();
    }
}

/**
 * @brief: Construct a basis function from a partition, index (global numbering) and curve.
 */
RegularP0_1D::BasisFunction::BasisFunction(const std::vector<double> &partition, int index, const Curve2D &curve) :
    P0Function(partition, BEM::basisVector(index, partition.size() - 1), curve)
{
    _support.reset(new BEM::Interval1D(_partition.at(index), _partition.at(index+1)));
}
