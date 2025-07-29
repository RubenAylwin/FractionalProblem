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
#include <Mesh.h>
#include <RegularP0Mesh.h>
#include <tbb/parallel_for.h>
#include <boost/math/quadrature/tanh_sinh.hpp>
useMessages("REG_MESH_P0");


/**
 * @brief: Constructor
 */
RegularP0Mesh_1D::RegularP0Mesh_1D(const Mesh1D &mesh) :
    DiscreteSpaceMesh(mesh),
    _size{_mesh.numElements()}
{
    _integrator1D.reset(new GaussLegendre_1D(30));
    _integrator2D.reset(new GaussLegendre_2D<5, 5>());
}

/**
 * @brief: Destructor
 */
RegularP0Mesh_1D::~RegularP0Mesh_1D(void)
{
}

/**
 * @brief: Test given function on the discretized curve against a specified element of the basis.
 */
BEM::Complex RegularP0Mesh_1D::testAgainstBaseElement(const ScalarFunctionBase_2D &function, unsigned globalNumber) const
{
    auto &basis = basisFunction(globalNumber);
    auto &element = _mesh.getElement(globalNumber);
    double jacobian = element.getSize();
    ExplicitScalarFunction_1D integrand([&basis, &function, &element, jacobian, globalNumber](double t) {
        return jacobian*basis.evaluate(globalNumber,t)*function(element(t).getX(), element(t).getY());
    });
    return _integrator1D->integrate(0, 1, integrand);
}


/**
 * @brief: Given a coefficient vector return the associated discrete function.
 * This is the only way to generate a discretization/construct a discrete function.
 */
std::unique_ptr<DiscreteFunctionMesh> RegularP0Mesh_1D::generateFunction(const std::vector<BEM::Complex> coefficients) const
{
    return std::unique_ptr<DiscreteFunctionMesh>(new P0Function(_mesh, coefficients));
}


/**
 * @brief: Given a function, return the discrete function with the same DoFs.
 */
std::unique_ptr<DiscreteFunctionMesh> RegularP0Mesh_1D::generateFunction(const ScalarFunctionBase_2D &function) const
{
    std::vector<BEM::Complex> coef(_size, 0.0);
    for (unsigned i = 0; i < _mesh.numElements(); i++) {
        ExplicitScalarFunction_1D integrand([&](double t){
            auto point = _mesh.getElement(i)(t);
            return function(point.getX(), point.getY());
        });
        coef[i] = _integrator1D->integrate(0, 1, integrand);
    }
    return generateFunction(coef);
}

/**
 * @brief: Returns the i-th basis function.
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const BasisFunctionMesh &RegularP0Mesh_1D::basisFunction(const unsigned globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber < _mesh.numElements() and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid basis index (out of bounds)");
    }

    std::lock_guard<std::mutex> lock(_basisMutex);
    if (_basis.find(globalNumber) == _basis.end()) {
        // Fill in coefficient vector
        std::vector<BEM::Complex> coefficients(_mesh.numElements(), std::complex(0.0, 0.0));

        // Coef = 1 for the corresponding basis
        coefficients[globalNumber] = std::complex(1.0, 0.0);

        // Return unique ptr
        _basis.insert({globalNumber ,std::unique_ptr<BasisFunctionMesh>(new BasisFunction(_mesh, globalNumber))});
    }
    return *(_basis.at(globalNumber));
}

/**
 * @brief: Construct a discrete function from given coefficients.
 */
RegularP0Mesh_1D::P0Function::P0Function(const Mesh1D &mesh, const std::vector<BEM::Complex> coefficients) :
    DiscreteFunctionMesh(mesh, coefficients)
{
    for (unsigned i = 0; i < _coefficients.size(); ++i) {
        auto &val = _coefficients[i];
        if (std::abs(val) > 1E-9) {
            _support.push_back(i);
        }
    }
    assert(_coefficients.size() == _mesh.numElements());
}

/**
 * @brief: Evaluate the discrete function on a given element and t specifying the point in the element
 * with a (0,1) discretization.
 */
BEM::Complex RegularP0Mesh_1D::P0Function::evaluate(unsigned elementIndex, [[maybe_unused]]double t) const
{
    assert(elementIndex < _mesh.numElements() and elementIndex >= 0);
    return _coefficients[elementIndex];
}

/**
 * @brief: Constructor of a Basis Function
 */
RegularP0Mesh_1D::BasisFunction::BasisFunction(const Mesh1D &mesh, unsigned index) :
    DiscreteFunctionMesh(mesh,BEM::basisVector(index, mesh.numElements())),
    BasisFunctionMesh(mesh,BEM::basisVector(index, mesh.numElements())),
    P0Function(mesh, BEM::basisVector(index, mesh.numElements()))
{
}
/**
 * @brief: return a function evaluating then derivative of the corresponding function. 
 */
DiscreteFunctionMesh &RegularP0Mesh_1D::P0Function::derivative(void) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (not _derivative) {
        std::vector<BEM::Complex> vec(_mesh.numElements(), 0.0);
        _derivative.reset(new P0Function(_mesh, vec));
    }
    return *_derivative;
}

/**
 * @brief: return a function evaluating the anti-derivative of the corresponding function. First corresponds
 * to the point of evaluation, second one to constant (since derivative of constant is 0).
 */
const ScalarFunctionBase_1D &RegularP0Mesh_1D::BasisFunction::antiDerivative(unsigned elementIndex) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (_antiDerivatives.find(elementIndex) == _antiDerivatives.end()) {
        auto antiDerivative = [] (double t) -> BEM::Complex {
            // Check for correct range (TODO: THIS SHOULD BE A THROW) rate
            if (t < 0 or t > 1) {
                msg(0) << "Got valid argument " << t << " not in [0, 1]." << endMsg;
                throw std::invalid_argument("Invalid argument, should be between 0 and 1");
            }
            return t;
        };
        _antiDerivatives.emplace(elementIndex, new ExplicitScalarFunction_1D(antiDerivative));
    }
    return *(_antiDerivatives[elementIndex]);
}
