#include <DiscreteSpace.h>
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
useMessages("DISC_SPACE");

/////////////////////////////////////////////////
// Right now this code is DEPRECATED in favour //
// of DiscreteSpaceMesh.                       //
/////////////////////////////////////////////////


/*
 * Code for Discrete space (base class).
 */

/**
 * @brief: Project a function onto the discrete space (L2 projection).
 */
BEM::ColVector DiscreteSpace_1D::projectFunction(const ScalarFunctionBase_1D &function) const
{
    auto matrix = identityOp();
    auto rhs = testAgainstBasis(function);

    BEM::ColVector solutionV = matrix.colPivHouseholderQr().solve(rhs);
    return solutionV;
}

/**
 * @brief: Given a coefficient vector, return the associated discrete function.
 * (override for call with a std::vector which is implemented on each discrete space.)
 */
DiscreteFunctionPtr DiscreteSpace_1D::generateFunction(const BEM::ColVector &coefficients) const
{
    std::vector<BEM::Complex> functionVector(coefficients.size(), 0);
    for (int i = 0; i < coefficients.size(); ++i) {
        functionVector[i] = coefficients[i];
    }
    return generateFunction(functionVector);
}

/**
 * @brief: Constructor.
 */
DiscreteSpace_1D::DiscreteSpace_1D(void) :
    _integrator1D(nullptr),
    _integrator2D(nullptr)
{
}

/**
 * @brief: Destructor.
 */
DiscreteSpace_1D::~DiscreteSpace_1D(void)
{
}

/**
 * @brief: Test a given function against a specified basis function (L2 product).
 */
BEM::Complex DiscreteSpace_1D::testAgainstBaseElement(const ScalarFunctionBase_1D &function, int globalNumber) const
{
    auto &basis = basisFunction(globalNumber);
    return _integrator1D->integrate(basis.brokenSupport(), *(basis*function));
}

/**
 * @brief: Given a scalar function, test it against the basis and return the vector.
 * @Assumption: The domain of the given function is [0, 1]
 */
BEM::ColVector DiscreteSpace_1D::testAgainstBasis(const ScalarFunctionBase_1D &function) const
{
    auto size = getSize();
    BEM::ColVector result(size);
    for (unsigned i = 0; i < size; ++i) {
        result[i] = testAgainstBaseElement(function, i);
    }
    return result;    
}

/**
 * @brief: Matrix of L2-products.
 */
BEM::Matrix DiscreteSpace_1D::identityOp(void) const
{
    auto size = getSize();
    BEM::Matrix matrix(size, size);
    for (unsigned i = 0; i < size; ++i) {
        auto &test = basisFunction(i);
        for (unsigned j = 0; j < size; ++j) {
            auto &trial = basisFunction(j);
            matrix(i, j) = _integrator1D->integrate(test.brokenSupport(), *(test*trial));
        }
    }
    return matrix;
}

/*
 * Code for P1 space on a curve
 */

/**
 * @brief: Constructor
 */
RegularP1_1D::RegularP1_1D(PolyPeriodicCurve &curve, int size, bool zeroBC) :
    RegularP1_1D(static_cast<Curve2D &>(curve), size, zeroBC)
{
    auto intrinsicPartition = curve.getPartition();
    assert(intrinsicPartition.size() < _size + 1U);

    for (const auto &val : intrinsicPartition) {
        if (val == 0 or val == 1) {
            continue;
        }
        double integerPart;
        double fractionalPart = std::modf(val*(_size-1), &integerPart);

        if (fractionalPart == 0) {
            continue;
        }
        int index = (fractionalPart > 0.5) ? integerPart + 1 : integerPart;
        _partition[index] = val;
    }
    
    // Correctness checks
    assert(_partition[0] == 0.0);
    assert(_partition[_size - 1] == 1.0);
}

/**
 * @brief: Constructor. Takes Curve, number of basis elements and a bool
 * (true if 0-boundary condition, false if no bc.)
 */
RegularP1_1D::RegularP1_1D(const Curve2D &curve, int size, bool zeroBC) :
    DiscreteSpaceOnCurve_1D(curve),
    _size(size),
    _zeroBC{zeroBC},
    _partition(_zeroBC ? _size + 2 : _size, 0.0),
    _baseName(_zeroBC ? "RegularP10" : "RegularP1")
{
    _integrator1D.reset(new GaussLegendre_1D(4));
    _integrator2D.reset(new GaussLegendre_2D<5, 6>());
    msg(0)<< "start::RegularP1_1D::FirstConstructor" << endMsg;
    // Size of each interval
    double intervalSize = 1.0/(_partition.size() - 1);
    msg(1) << "Filling nodes" << endMsg;
    // Fill in the nodes of the partition
    for (size_t i = 0; i < _partition.size(); ++i) {
        _partition[i] = i*intervalSize;
    }
    _elements.emplace_back(0, 1, 0, _zeroBC ? 1 : 2);
    size_t i = 1;
    for (; i < _partition.size() - 2; ++i) {
        _elements.emplace_back(i, i+1, i, 2);
    }
    _elements.emplace_back(i, i+1, i, _zeroBC ? 1 : 2);
    msg(1) << "Begin checks" << endMsg;
    // Correctness checks
    assert(_partition[0] == 0.0);
    // assert(_partition[_size - 1] == 1);
    msg(0)<< "end::RegularP1_1D::FirstConstructor" << endMsg;
}

/**
 * @brief: Destructor
 */
RegularP1_1D::~RegularP1_1D(void)
{
}

/**
 * @brief: Returns the i-th basis function.
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const BasisFunction_1D &RegularP1_1D::basisFunction(const int globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber < static_cast<int>(_size) and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid basis index (out of bounds)");
    }

    std::lock_guard<std::mutex> lock(_basisMutex);
    if (_basis.find(globalNumber) == _basis.end()) {
        int localNumber = _zeroBC ? globalNumber + 1 : globalNumber;

        // Fill in coefficient vector
        std::vector<BEM::Complex> coefficients(_zeroBC ? _size + 2 : _size, std::complex(0.0, 0.0));

        // Coef = 1 for the corresponding basis
        coefficients[localNumber] = std::complex(1.0, 0.0);

        // Return unique ptr
        _basis.insert({globalNumber ,std::unique_ptr<BasisFunction_1D>(new BasisFunction(_partition, localNumber, _curve))});
    }
    return *(_basis.at(globalNumber));
}

/**
 * @brief: Returns the i-th DOF.
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const DiscreteFunctional_1D &RegularP1_1D::degreeOfFreedom(const int globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber < static_cast<int>(_size) and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid DoF index (out of bounds)");
    }
    
    std::lock_guard<std::mutex> lock(_dofMutex);
    int localNumber = globalNumber;
    if (_zeroBC) {
        localNumber += 1;
    }
    if (_dof.find(globalNumber) == _dof.end()) {
        _dof.insert({globalNumber ,std::unique_ptr<DiscreteFunctional_1D>(new DiscreteFunctional_1D([this, localNumber](const ScalarFunctionBase_1D &function) -> BEM::Complex {
            return function(_partition[localNumber])*_curve.jacobian(_partition[localNumber]);
        }))});
    }
    return *(_dof.at(globalNumber));
}

/**
 * @brief: Given an element and local function numbering, return the global numbering.
 */
unsigned RegularP1_1D::globalIndex(unsigned element, unsigned localFunction) const
{
    assert(localFunction == 0 or localFunction == 1);
    assert(element < _zeroBC ? _size + 2 : _size);
    if (element == 0) {
        return localFunction;
    }
    
    return element + localFunction - (_zeroBC ? 1 : 0);
}

/**
 * @brief: Given a coefficient vector return the associated discrete function.
 * This is the only way to generate a discretization/construct a discrete function.
 */
std::unique_ptr<DiscreteFunction_1D> RegularP1_1D::generateFunction(const std::vector<BEM::Complex> coefficients) const
{
    msg(1) << "start::RegularP1_1D::generateFunction" << endMsg;
    assert(static_cast<unsigned>(coefficients.size()) == _size);
    msg(1) << "end::RegularP1_1D::generateFunction" << endMsg;
    return std::unique_ptr<DiscreteFunction_1D>(new P1Function(_partition, coefficients, _curve));
}

/**
 * @brief: Construct a discrete function from given coefficients.
 */
RegularP1_1D::P1Function::P1Function(const std::vector<double> &partition, const std::vector<BEM::Complex> coefficients, const Curve2D &curve) :
    _partition(partition),
    _coefficients(coefficients),
    _curve(curve)
{
    msg(1) << "start::P1Function::Constructor" << endMsg;
    if (_coefficients.size() == _partition.size() - 2) {
        _coefficients.insert(_coefficients.begin(), BEM::Complex(0.0));
        _coefficients.insert(_coefficients.end(), BEM::Complex(0.0));
    }
    msg(1) << "P1Function::Constructor Check" << endMsg;
    assert(_coefficients.size() == _partition.size());
    msg(1) << "end::P1Function::Constructor" << endMsg;
}

/**
 * @brief: Evaluate the function.
 */
BEM::Complex RegularP1_1D::P1Function::operator()(const double t) const
{
    // Check for correct range (TODO: THIS SHOULD BE A THROW)
    assert(t >= 0 and t <= 1);

    // Find the index
    auto index = BEM::findIndex(_partition, t, BasisFunction_1D::Direction::LEFT);

    // Get the coeffs
    auto leftCoeff = _coefficients[index];
    auto rightCoeff = _coefficients[index+1];

    auto leftPart = _partition[index];
    auto rightPart = _partition[index+1];
    // Return
    return (leftCoeff + (t - leftPart)*(rightCoeff-leftCoeff)/(rightPart-leftPart))/_curve.jacobian(t);
}


/**
 * @brief: Return the support of a given function
 */
BEM::Interval1D RegularP1_1D::P1Function::support(void) const
{
    if (_support) {
        return *_support;
    }

    // Find first non zero coefficient
    unsigned i = 0u;
    while (i < _coefficients.size() and std::abs(_coefficients[i]) < 1E-10) {
        ++i;
    }

    // If got to the end of the list, support is null.
    if (i == _coefficients.size()) {
        return {0,0};
    }

    // Correctness check
    assert(i < _coefficients.size());
    unsigned lower = i > 0 ? i - 1 : 0;
    
    // Find the last non zero coefficient
    i = _coefficients.size() - 1;
    while (i > lower and std::abs(_coefficients[i]) < 1E-10) {
        --i;
    }
    unsigned upper = i < _coefficients.size() - 1 ? i + 1 : i;

    
    // return values
    _support.reset(new BEM::Interval1D(_partition.at(lower), _partition.at(upper)));

    _brokenSupport.reset(new BEM::Support1DL{BEM::Interval1D(_partition.at(lower), _partition.at(lower + 1))});
    lower++;
    while (lower < upper) {
        _brokenSupport->push_back(BEM::Interval1D(_partition.at(lower), _partition.at(lower + 1)));
        lower++;
    }
    
    return *_support;
}

/**
 * @brief: Return the support of a given function
 */
const BEM::Support1DL &RegularP1_1D::P1Function::brokenSupport(void) const
{
    if (_brokenSupport) {
        return *_brokenSupport;
    }
    support();
    assert(_brokenSupport);
    return *_brokenSupport;
}

/**
 * @brief: return a function evaluating then derivative of the corresponding function. 
 */
ScalarFunctionBase_1D *RegularP1_1D::BasisFunction::derivative(void) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (not _derivative) {
        auto supp = brokenSupport();
        assert(supp.size() >= 1 and supp.size() <= 2);
        double a = std::min(supp.front().first, supp.back().first);
        double c = std::max(supp.front().second, supp.back().second);
        double b = supp.size() == 2 ? std::max(supp.front().first, supp.back().first) : (a + c)/2;
        auto va = this->operator()(a);
        auto vb = this->operator()(b);
        auto vc = this->operator()(c);
        auto derivative = [a, b, c, va, vb, vc](double t) -> BEM::Complex {
            if (t < a or t > c) {
                return 0.0;
            }

            if (t < b) {
                return (vb - va)/(b-a);
            }

            return (vc - vb)/(c-b);    
        };
            _derivative.reset(new ExplicitScalarFunction_1D(derivative));
    }
    
    return _derivative.get();
}

/**
 * @brief: return a function evaluating the anti-derivative of the corresponding function. First corresponds
 * to the point of evaluation, second one to constant (since derivative of constant is 0).
 */
ScalarFunctionBase_2D *RegularP1_1D::BasisFunction::antiDerivative(Direction direction) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    
    if (direction == Direction::LEFT) {
        if (not _antiDerivativeLeft) {
            if (_index==0) {
                auto antiDerivativeLeft = [&] (double t, double s) -> BEM::Complex {
                    return (t+s) - (t+s)*(t-s)*0.5/_partition[1];
                };
                _antiDerivativeLeft.reset(new ExplicitScalarFunction_2D(antiDerivativeLeft));
            } else if (_index == _partition.size() - 1U ) {
                auto antiDerivativeLeft = [&] (double t, double s) -> BEM::Complex {
                    return ((t+s)*(t-s)*0.5 - (t+s)*_partition[_index-1])/(1.-_partition[_index-1]);
                };
                _antiDerivativeLeft.reset(new ExplicitScalarFunction_2D(antiDerivativeLeft));
            } else  {
                auto antiDerivativeLeft = [&] (double t, double s) -> BEM::Complex {
                    double a = _partition[_index-1];
                    double b = _partition[_index];
                    double c = _partition[_index+1];

                    double value = 0;
                    if (t <= b) {
                        value = (0.5*t*t - a*t)/(b-a);
                    } else {
                        double valbL = (0.5*b*b - a*b)/(b-a);
                        double valbR = (-0.5*b*b + c*b)/(c-b);
                        value = (-0.5*t*t + c*t)/(c-b) + valbL - valbR;
                    }
                    double constant = 0;
                    if (-s <= b) {
                        constant = (-0.5*s*s - a*s)/(b-a);
                    } else {
                        double valbL = (0.5*b*b - a*b)/(b-a);
                        double valbR = (-0.5*b*b + c*b)/(c-b);
                        constant = (0.5*s*s + c*s)/(c-b) - valbL + valbR;
                    }
                    return value + constant;
                };
                _antiDerivativeLeft.reset(new ExplicitScalarFunction_2D(antiDerivativeLeft));
            }
        }
        return _antiDerivativeLeft.get();

    } else /*if (direction == Direction::RIGHT)*/ {
        if (not _antiDerivativeRight) {
            if (_index==0) {
                auto antiDerivativeRight = [&] (double t, double s) -> BEM::Complex {
                    return (t+s) - (t+s)*(t-s)*0.5/_partition[1];
                };
                _antiDerivativeRight.reset(new ExplicitScalarFunction_2D(antiDerivativeRight));
            } else if (_index==_partition.size()-1) {
                auto antiDerivativeRight = [&] (double t, double s) -> BEM::Complex {
                    return ((t+s)*(t-s)*0.5 - (t+s)*_partition[_index-1])/(1.-_partition[_index-1]);
                };
                _antiDerivativeRight.reset(new ExplicitScalarFunction_2D(antiDerivativeRight));
            } else  {
                auto antiDerivativeRight = [&] (double t, double s) -> BEM::Complex {
                    double a = _partition[_index-1];
                    double b = _partition[_index];
                    double c = _partition[_index+1];

                    double value = 0;
                    if (t < b) {
                        value = (0.5*t*t - a*t)/(b-a);
                    } else {
                        double valbL = (0.5*b*b - a*b)/(b-a);
                        double valbR = (-0.5*b*b + c*b)/(c-b);
                        value = (-0.5*t*t + c*t)/(c-b) + valbL - valbR;
                    }
                    double constant = 0;
                    if (-s < b) {
                        constant = (-0.5*s*s - a*s)/(b-a);
                    } else {
                        double valbL = (0.5*b*b - a*b)/(b-a);
                        double valbR = (-0.5*b*b + c*b)/(c-b);
                        constant = (0.5*s*s + c*s)/(c-b) - valbL + valbR;
                    }
                    return value + constant;
                };
                _antiDerivativeRight.reset(new ExplicitScalarFunction_2D(antiDerivativeRight));
            }
        }
        return _antiDerivativeRight.get();
    }
}


/**
 * @brief: Construction of a basis function on a given partition and curve for the specified index.
 */
RegularP1_1D::BasisFunction::BasisFunction(const std::vector<double> &partition, unsigned index, const Curve2D &curve) :
    P1Function(partition, BEM::basisVector(index, partition.size()), curve),
    _index{index}
{
    if (index == 0) {
        _support.reset(new BEM::Interval1D(_partition.at(index), _partition.at(index+1)));
        _brokenSupport.reset(new BEM::Support1DL{*_support});
    } else if (index == _partition.size() - 1) {
        _support.reset(new BEM::Interval1D(_partition.at(index-1), _partition.at(index)));
        _brokenSupport.reset(new BEM::Support1DL{*_support});
    } else {
        _support.reset(new BEM::Interval1D(_partition.at(index-1), _partition.at(index+1)));
        _brokenSupport.reset(new BEM::Support1DL{BEM::Interval1D(_partition.at(index-1), _partition.at(index)), BEM::Interval1D(_partition.at(index), _partition.at(index+1))});
    }
}

/**
 * @brief: Constructor on a periodic curve.
 */
GeoP1_1D::GeoP1_1D(PolyPeriodicCurve &curve, int size, double power, bool zeroBC) :
    RegularP1_1D(curve, size, zeroBC),
    _power{power}
{
}

/**
 * @brief: Constructor on a general curve.
 */
GeoP1_1D::GeoP1_1D(const Curve2D &curve, int size, double power, bool zeroBC) :
    RegularP1_1D(curve, size, zeroBC),
    _power{power}
{
    for (auto &val : _partition) {
        val = std::pow(val, _power);
    }
}
