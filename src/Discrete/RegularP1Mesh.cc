#include <unordered_set>
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
#include <RegularP1Mesh.h>
#include <RegularP0Mesh.h>

useMessages("REG_MESH_P1");
/**
 * @brief: Constructor
 */
RegularP1_0Mesh_1D::RegularP1_0Mesh_1D(const Mesh1D &mesh) :
    DiscreteSpaceMesh(mesh),
    _size{_mesh.numElements() - 1}
{
    _integrator1D.reset(new GaussLegendre_1D(2));
    _integrator2D.reset(new GaussLegendre_2D<5, 5>());
}

/**
 * @brief: Destructor
 */
RegularP1_0Mesh_1D::~RegularP1_0Mesh_1D(void)
{
}

/**
 * @brief: Testing against a base element.
 */
BEM::Complex RegularP1_0Mesh_1D::testAgainstBaseElement(const ScalarFunctionBase_2D &function, unsigned globalNumber) const
{
    // For this space, all basis functions are supported on two elements
    auto &basis = basisFunction(globalNumber);
    auto &element1 = _mesh.getElement(globalNumber);
    auto &element2 = _mesh.getElement(globalNumber+1);
    double jacobian1 = element1.getSize();
    double jacobian2 = element2.getSize();
    auto integrand = [&basis, &function, &element1, &element2, jacobian1, jacobian2, globalNumber](double t) {
        return jacobian1*basis.evaluate(globalNumber,t)*function(element1(t).getX(), element1(t).getY()) +
            jacobian2*basis.evaluate(globalNumber+1,t)*function(element2(t).getX(), element2(t).getY());
    };
    boost::math::quadrature::tanh_sinh<double> integratorTS;
    return integratorTS.integrate(integrand, 0, 1);
    ExplicitScalarFunction_1D integrandE(std::move(integrand));
    return _integrator1D->integrate(0, 1, integrandE);    
}

/**
 * @brief: Given a coefficient vector return the associated discrete function.
 * This is the only way to generate a discretization/construct a discrete function.
 */
std::unique_ptr<DiscreteFunctionMesh> RegularP1_0Mesh_1D::generateFunction(const std::vector<BEM::Complex> coefficients) const
{
    return std::unique_ptr<DiscreteFunctionMesh>(new P1Function(_mesh, coefficients, *_integrator1D));
}

/**
 * @brief: Given a function, return the function with the same dof
 */
std::unique_ptr<DiscreteFunctionMesh> RegularP1_0Mesh_1D::generateFunction(const ScalarFunctionBase_2D &function) const
{
    std::vector<BEM::Complex> coef(_size, 0.0);
    for (unsigned i = 1; i < _mesh.numElements(); i++) {
        auto point = _mesh.point(i);
        coef[i-1] = function(point.getX(), point.getY());
        
    }
    return generateFunction(coef);
}

/**
 * @brief: Returns the i-th basis function.
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const BasisFunctionMesh &RegularP1_0Mesh_1D::basisFunction(const unsigned globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber < _mesh.numElements() - 1 and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid basis index (out of bounds)");
    }

    std::lock_guard<std::mutex> lock(_basisMutex);
    if (_basis.find(globalNumber) == _basis.end()) {
        // Fill in coefficient vector
        std::vector<BEM::Complex> coefficients(_size, std::complex(0.0, 0.0));

        // Coef = 1 for the corresponding basis
        coefficients[globalNumber] = std::complex(1.0, 0.0);

        // Return unique ptr
        _basis.insert({globalNumber ,std::unique_ptr<BasisFunctionMesh>(new BasisFunction(_mesh, globalNumber, *_integrator1D))});
    }
    return *(_basis.at(globalNumber));
}

/**
 * @brief: Construct a discrete function from given coefficients.
 */
RegularP1_0Mesh_1D::P1Function::P1Function(const Mesh1D &mesh, const std::vector<BEM::Complex> coefficients, Integrator_1D &integrator) :
    DiscreteFunctionMesh(mesh, coefficients),
    _size{_mesh.numElements() - 1},
    _integrator{integrator}
{
    assert(_coefficients.size() == _size);
    std::unordered_set<unsigned> support;
    for (unsigned i = 0; i < _coefficients.size(); ++i) {
        auto &val = _coefficients[i];
        if (std::abs(val) > 1E-9) {
            support.insert(i);
            support.insert(i+1);
        }
    }
    for (auto &s : support) {
        _support.push_back(s);
    }
    
}

/**
 * @brief: Evaluate the function on a specified element and t a point in the element assuming a linear
 * [0,1] parametrization.
 */
BEM::Complex RegularP1_0Mesh_1D::P1Function::evaluate(unsigned elementIndex, [[maybe_unused]]double t) const
{
    assert(elementIndex < _mesh.numElements() and elementIndex >= 0);
    BEM::Complex a = elementIndex > 0 ? _coefficients[elementIndex - 1] : 0.0;
    BEM::Complex b = elementIndex < _size ? _coefficients[elementIndex] : 0.0;
    return a+(b-a)*t;
}

/**
 * @brief: return a function evaluating then derivative of the corresponding function. 
 */
DiscreteFunctionMesh &RegularP1_0Mesh_1D::P1Function::derivative(void) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (not _derivative) {
        auto vec1 = _coefficients;
        auto vec2 = _coefficients;
        vec1.push_back(0.0);
        vec2.insert(vec2.begin(), 0.0);
        auto coeff = (vec1 - vec2);
        for (unsigned i = 0; i < coeff.size(); i++) {
            coeff[i]/=_mesh.getElement(i).getSize();
        }
        _derivative.reset(new RegularP0Mesh_1D::P0Function(_mesh,coeff));
    }
    return *_derivative;
}

/**
 * @brief: Return a function evaluating the RL derivative of the given function.
 * @input: Order is a number between -100 and 100 specifying the order of the derivative*100.
 * If positive->left derivative, if negative->right derivative.
 * @output: Reference to a discrete P0 function on the mesh that is the projection of the derivative.
 * @desc: RL derivative is computed using FT of calculus. The RL derivative is (d/dx)RL_Int. The P0
 * projection depends on the integral over the element. Therefore, only on the evaluation of the RL_Int
 * on the limits of the element.
 */
DiscreteFunctionMesh &RegularP1_0Mesh_1D::P1Function::derivative(int order) const
{
    // Check that the order is in the correct range.
    assert(order < 100 and order > -100 and order != 0);

    // Correct order
    double fOrder = std::abs(order*0.01);

    //Expect a parallell call, so lock the function.
    std::lock_guard<std::mutex> lock(_derMutex);

    //If the derivative has not been constructed beforehand, we have to build it.
    if (_fracDerivatives.find(order) == _fracDerivatives.end()) {
        // Construct the TS integrator only if necessary.
        if ((not _integratorTS.has_value())) {
            _integratorTS.emplace(boost::math::quadrature::tanh_sinh<double>());
        }
        
        //coefficient vector for the 
        std::vector<BEM::Complex> coef(_mesh.numElements(), 0.0);
        // if order is positive, left derivative
        if (order > 0) {
            // iterate over elements
            for (unsigned i = 0; i < _mesh.numElements(); ++i) {
                if (i < _support[0]) {
                    coef[i] = 0.0; //zero initialization.
                }
                auto evaluationElement = _mesh.getElement(i);
                //Integration over the same element.
                auto integrandSelf = [&](double t, double xc) -> double {
                    if (t > 0.9) {
                        return evaluate(i, t).real()*std::pow(xc*evaluationElement.getSize(), -fOrder);
                    }
                    return evaluate(i, t).real()*std::pow(evaluationElement(1.).getX()-evaluationElement(t).getX(), -fOrder);
                };
                coef[i] = _integratorTS->integrate(integrandSelf, 0, 1)*evaluationElement.getSize();
                // Integration over all previous elements
                for (unsigned j = 0; j < i; ++j) {
                    auto integrationElement = _mesh.getElement(j);
                    // integrand against point B
                    auto integrandB = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((evaluationElement(1.) - integrationElement(t)).norm(), -fOrder);
                    };
                    // integrand against point A
                    auto integrandA = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((evaluationElement(0.) - integrationElement(t)).norm(), -fOrder);
                    };
                    // If elements are far away, integrand against point A is not singular.
                    if (j < i - 1) {
                        coef[i] += (_integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandB))) - _integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandA))))*integrationElement.getSize();
                    } else {
                        // If the elements are close, the integrand against point A is singular, use TS integration.
                        auto integrandA2 = [&](double t, double xc) -> double {
                            if (t > 0.9) {
                                return evaluate(j, 1.-xc).real()*std::pow(xc*integrationElement.getSize(), -fOrder);
                            }
                            return evaluate(j, t).real()*std::pow((evaluationElement(0.) - integrationElement(t)).norm(), -fOrder);
                        };
                        coef[i] += (_integrator.integrate(0, 1., ExplicitScalarFunction_1D(std::move(integrandB))) - _integratorTS->integrate(integrandA2, 0., 1.))*integrationElement.getSize();
                    }
                }
                // Remaining scale factors
                coef[i] = coef[i]/(std::tgamma(1.-fOrder)*evaluationElement.getSize());
            }
        } else {
            // Right RL derivative. Analogous as left.
            for (unsigned i = 0; i < _mesh.numElements(); ++i) {
                if (i > _support.back()) {
                    coef[i] = 0.0;
                }
                auto evaluationElement = _mesh.getElement(i);
                auto integrandSelf = [&](double t, double xc) -> double {
                    if (t < 0.1) {
                        return evaluate(i, t).real()*std::pow(-xc*evaluationElement.getSize(), -fOrder);
                    }
                    return evaluate(i, t).real()*std::pow(evaluationElement(t).getX()-evaluationElement(0).getX(), -fOrder);
                };
                coef[i] = _integratorTS->integrate(integrandSelf, 0, 1)*evaluationElement.getSize();
                int int_i(i);
                for (int j = _mesh.numElements() - 1; j > int_i; --j) {
                    auto integrationElement = _mesh.getElement(j);
                    auto integrandA = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((integrationElement(t) - evaluationElement(0.)).norm(), -fOrder);
                    };
                    auto integrandB = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((integrationElement(t) - evaluationElement(1.)).norm(), -fOrder);
                    };

                    if (j > int_i + 1) {
                        coef[i] -= (_integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandB))) - _integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandA))))*integrationElement.getSize();
                    } else {
                        auto integrandB = [&](double t, double xc) -> double {
                            if (t < 0.1) {
                                return evaluate(j, t).real()*std::pow(-xc*integrationElement.getSize(), -fOrder);
                        }
                            return evaluate(j, t).real()*std::pow((integrationElement(t) - evaluationElement(1.)).norm(), -fOrder);
                    };
                        coef[i] -= (_integratorTS->integrate(integrandB, 0., 1.) - _integrator.integrate(0, 1., ExplicitScalarFunction_1D(std::move(integrandA))))*integrationElement.getSize();
                    }
                }
                coef[i] = coef[i]/(std::tgamma(1.-fOrder)*evaluationElement.getSize());
            }            
        }
        _fracDerivatives.emplace(order, new RegularP0Mesh_1D::P0Function(_mesh,coef));
    }
    return *(_fracDerivatives[order]);
}

/**
 * @brief: Construct a basis function
 */
RegularP1_0Mesh_1D::BasisFunction::BasisFunction(const Mesh1D &mesh, unsigned index, Integrator_1D & integrator) :
    DiscreteFunctionMesh(mesh, BEM::basisVector(index, mesh.numElements() - 1)),
    BasisFunctionMesh(mesh, BEM::basisVector(index, mesh.numElements() - 1)),
    P1Function(mesh, BEM::basisVector(index, mesh.numElements() - 1), integrator)
{
}

/**
 * @brief: return a function evaluating the anti-derivative of the corresponding function. First corresponds
 * to the point of evaluation, second one to constant (since derivative of constant is 0).
 */
const ScalarFunctionBase_1D &RegularP1_0Mesh_1D::BasisFunction::antiDerivative(unsigned elementIndex) const
{
    assert(false and "Not Implemented");
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

/**
 * @brief: Constructor
 */
RegularP1_Mesh_1D::RegularP1_Mesh_1D(const Mesh1D &mesh) :
    DiscreteSpaceMesh(mesh),
    _size{_mesh.numElements() + 1}
{
    _integrator1D.reset(new GaussLegendre_1D(2));
    _integrator2D.reset(new GaussLegendre_2D<5, 5>());
}

/**
 * @brief: Destructor
 */
RegularP1_Mesh_1D::~RegularP1_Mesh_1D(void)
{
}

BEM::Complex RegularP1_Mesh_1D::testAgainstBaseElement(const ScalarFunctionBase_2D &function, unsigned globalNumber) const
{
    auto &basis = basisFunction(globalNumber);
    BEM::Complex result = 0.0;
    for (auto &i : basis.support()) {
        auto &element = _mesh.getElement(i);
        double jacobian = element.getSize();
        auto integrand = [&basis, &function, &element, jacobian, i](double t) {
            return jacobian*basis.evaluate(i,t)*function(element(t).getX(), element(t).getY());
        };
        boost::math::quadrature::tanh_sinh<double> integratorTS;
        result += integratorTS.integrate(integrand, 0, 1);
    }
    return result;
}


/**
 * @brief: Given a coefficient vector return the associated discrete function.
 * This is the only way to generate a discretization/construct a discrete function.
 */
std::unique_ptr<DiscreteFunctionMesh> RegularP1_Mesh_1D::generateFunction(const std::vector<BEM::Complex> coefficients) const
{
    return std::unique_ptr<DiscreteFunctionMesh>(new P1Function(_mesh, coefficients, *_integrator1D));
}

/**
 * @brief: Given a function, return the function with the same dof
 */
std::unique_ptr<DiscreteFunctionMesh> RegularP1_Mesh_1D::generateFunction(const ScalarFunctionBase_2D &function) const
{
    std::vector<BEM::Complex> coef(_size, 0.0);
    for (unsigned i = 0; i < _size; i++) {
        auto point = _mesh.point(i);
        coef[i] = function(point.getX(), point.getY());
    }
    return generateFunction(coef);
}

/**
 * @brief: Returns the i-th basis function.
 */
/*TODO: CONSIDER RETURNING POINTER TO SCALARFUNCTIONBASE_1D*/
const BasisFunctionMesh &RegularP1_Mesh_1D::basisFunction(const unsigned globalNumber) const
{
    // Check numbering. Throw exception if out of bounds.
    if (not (globalNumber <= _mesh.numElements() and globalNumber >= 0)) {
        throw std::invalid_argument("Invalid basis index (out of bounds)");
    }

    std::lock_guard<std::mutex> lock(_basisMutex);
    if (_basis.find(globalNumber) == _basis.end()) {
        // Fill in coefficient vector
        std::vector<BEM::Complex> coefficients(_size, std::complex(0.0, 0.0));

        // Coef = 1 for the corresponding basis
        coefficients[globalNumber] = std::complex(1.0, 0.0);

        // Return unique ptr
        _basis.insert({globalNumber ,std::unique_ptr<BasisFunctionMesh>(new BasisFunction(_mesh, globalNumber, *_integrator1D))});
    }
    return *(_basis.at(globalNumber));
}


/*
 * Code for the discrete P0 functions.
 */

/**
 * @brief: Construct a discrete function from given coefficients.
 */
RegularP1_Mesh_1D::P1Function::P1Function(const Mesh1D &mesh, const std::vector<BEM::Complex> coefficients, Integrator_1D &integrator) :
    DiscreteFunctionMesh(mesh, coefficients),
    _size{_mesh.numElements() + 1},
    _integrator{integrator}
{
    assert(_coefficients.size() == _size);
    std::unordered_set<unsigned> support;
    for (unsigned i = 0; i < _size; ++i) {
        auto &val = _coefficients[i];
        if (std::abs(val) > 1E-9) {
            if (i > 0 and i < _mesh.numElements()) {
                support.insert(i);
                support.insert(i-1);
            } else if (i==0) {
                support.insert(i);
            } else {
                support.insert(i-1);
            }
        }
    }
    for (auto &s : support) {
        _support.push_back(s);
    }
    
}

/**
 * @brief: Evaluate the function on a specified element and t a point in the element assuming a linear
 * [0,1] parametrization.
 */
BEM::Complex RegularP1_Mesh_1D::P1Function::evaluate(unsigned elementIndex, [[maybe_unused]]double t) const
{
    assert(elementIndex < _mesh.numElements() +1 and elementIndex >= 0);
    BEM::Complex a = _coefficients[elementIndex];
    BEM::Complex b = _coefficients[elementIndex+1];
    return a+(b-a)*t;
}

/**
 * @brief: return a function evaluating then derivative of the corresponding function. 
 */
DiscreteFunctionMesh &RegularP1_Mesh_1D::P1Function::derivative(void) const
{
    std::lock_guard<std::mutex> lock(_derMutex);
    if (not _derivative) {
        std::vector<BEM::Complex> coeff(_mesh.numElements(), 0.0);
        for (unsigned i = 0; i < _mesh.numElements(); i++) {
            coeff[i]=(_coefficients[i+1]-_coefficients[i])/_mesh.getElement(i).getSize();
        }
        _derivative.reset(new RegularP0Mesh_1D::P0Function(_mesh,coeff));
    }
    return *_derivative;
}


/**
 * @brief: return a function evaluating the caputo derivative of the function.
 * @desc: Computed by computing the RL derivative and subtracting the singular term.
 * RL computation as before plus a term taking care of the singular part.
 */
DiscreteFunctionMesh &RegularP1_Mesh_1D::P1Function::derivative(int order) const
{
    assert(order < 100 and order > -100 and order != 0);
    double fOrder = std::abs(order*0.01);
    std::lock_guard<std::mutex> lock(_derMutex);
    if (_fracDerivatives.find(order) == _fracDerivatives.end()) {
        if ((not _integratorTS.has_value())) {
            _integratorTS.emplace(boost::math::quadrature::tanh_sinh<double>());
        }
        std::vector<BEM::Complex> coef(_mesh.numElements(), 0.0);
        if (order > 0) {
            for (unsigned i = 0; i < _mesh.numElements(); ++i) {
                if (i < _support[0]) {
                    coef[i] = 0.0;
                }
                auto evaluationElement = _mesh.getElement(i);
                auto integrandSelf = [&](double t, double xc) -> double {
                    if (t > 0.9) {
                        return evaluate(i, t).real()*std::pow(xc*evaluationElement.getSize(), -fOrder);
                    }
                    return evaluate(i, t).real()*std::pow(evaluationElement(1.).getX()-evaluationElement(t).getX(), -fOrder);
                };
                auto xA = evaluationElement.getA().getX();
                auto xB = evaluationElement.getB().getX();
                coef[i] = _integratorTS->integrate(integrandSelf, 0, 1)*evaluationElement.getSize() - evaluate(0, 0)*(std::pow(xB, 1.-fOrder)-std::pow(xA, 1.-fOrder))/(std::tgamma(1.-fOrder)*(1.-fOrder));
                int int_i(i);
                for (int j = 0; j < int_i; ++j) {
                    auto integrationElement = _mesh.getElement(j);
                    auto integrandA = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((evaluationElement(0.) - integrationElement(t)).norm(), -fOrder);
                    };
                    auto integrandB = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((evaluationElement(1.) - integrationElement(t)).norm(), -fOrder);
                    };
                    if (j < int_i - 1) {
                        coef[i] += (_integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandB))) - _integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandA))))*integrationElement.getSize();
                    } else {
                        auto integrandA2 = [&](double t, double xc) -> double {
                            if (t > 0.9) {
                                return evaluate(j, 1.-xc).real()*std::pow(xc*integrationElement.getSize(), -fOrder);
                            }
                            return evaluate(j, t).real()*std::pow((evaluationElement(0.) - integrationElement(t)).norm(), -fOrder);
                        };
                        coef[i] += (_integrator.integrate(0, 1., ExplicitScalarFunction_1D(std::move(integrandB))) - _integratorTS->integrate(integrandA2, 0., 1.))*integrationElement.getSize();
                    }
                }
                coef[i] = coef[i]/(std::tgamma(1.-fOrder)*evaluationElement.getSize());

            }
        } else {
            for (unsigned i = 0; i < _mesh.numElements(); ++i) {
                if (i > _support.back()) {
                    coef[i] = 0.0;
                }
                auto evaluationElement = _mesh.getElement(i);
                auto integrandSelf = [&](double t, double xc) -> double {
                    if (t < 0.1) {
                        return evaluate(i, t).real()*std::pow(-xc*evaluationElement.getSize(), -fOrder);
                    }
                    return evaluate(i, t).real()*std::pow(evaluationElement(t).getX()-evaluationElement(0).getX(), -fOrder);
                };
                auto xA = evaluationElement.getA().getX();
                auto xB = evaluationElement.getB().getX();
                coef[i] = _integratorTS->integrate(integrandSelf, 0, 1)*evaluationElement.getSize() - evaluate(_mesh.numElements()-1., 1.)*(std::pow(1.-xA, 1.-fOrder)-std::pow(std::abs(1.-xB), 1.-fOrder))/(std::tgamma(1.-fOrder)*(1.-fOrder));
                int int_i(i);
                for (int j = _mesh.numElements() - 1; j > int_i; --j) {
                    auto integrationElement = _mesh.getElement(j);
                    auto integrandA = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((integrationElement(t) - evaluationElement(0.)).norm(), -fOrder);
                    };
                    auto integrandB = [&](double t) -> double {
                        return evaluate(j, t).real()*std::pow((integrationElement(t) - evaluationElement(1.)).norm(), -fOrder);
                    };

                    if (j > int_i + 1) {
                        coef[i] -= (_integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandB))) - _integrator.integrate(0., 1., ExplicitScalarFunction_1D(std::move(integrandA))))*integrationElement.getSize();
                    } else {
                        auto integrandB = [&](double t, double xc) -> double {
                            if (t < 0.1) {
                                return evaluate(j, t).real()*std::pow(-xc*integrationElement.getSize(), -fOrder);
                        }
                            return evaluate(j, t).real()*std::pow((integrationElement(t) - evaluationElement(1.)).norm(), -fOrder);
                    };
                        coef[i] -= (_integratorTS->integrate(integrandB, 0., 1.) - _integrator.integrate(0, 1., ExplicitScalarFunction_1D(std::move(integrandA))))*integrationElement.getSize();
                    }
                }
                coef[i] = coef[i]/(std::tgamma(1.-fOrder)*evaluationElement.getSize());
            }            
        }
        _fracDerivatives.emplace(order, new RegularP0Mesh_1D::P0Function(_mesh,coef));
    }
    return *(_fracDerivatives[order]);
}

/**
 * @brief: Constructor
 */
RegularP1_Mesh_1D::BasisFunction::BasisFunction(const Mesh1D &mesh, unsigned index, Integrator_1D & integrator) :
    DiscreteFunctionMesh(mesh, BEM::basisVector(index, mesh.numElements() + 1)),
    BasisFunctionMesh(mesh, BEM::basisVector(index, mesh.numElements() + 1)),
    P1Function(mesh, BEM::basisVector(index, mesh.numElements() + 1), integrator)
{
}

/**
 * @brief: return a function evaluating the anti-derivative of the corresponding function. First corresponds
 * to the point of evaluation, second one to constant (since derivative of constant is 0).
 */
const ScalarFunctionBase_1D &RegularP1_Mesh_1D::BasisFunction::antiDerivative(unsigned elementIndex) const
{
    assert(false and "Not Implemented");
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
