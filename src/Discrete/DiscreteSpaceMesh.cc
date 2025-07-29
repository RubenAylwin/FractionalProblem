#include <DiscreteSpaceMesh.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <ScalarValuedFunction.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <MyTypes.h>
#include <Utilities.h>
#include <Msg.h>
#include <boost/math/quadrature/tanh_sinh.hpp>
useMessages("DISC_SPACE_MESH");


/**
 * @brief: Constructor.
 */
DiscreteSpaceMesh::DiscreteSpaceMesh(const Mesh1D &mesh) :
    _mesh{mesh},
    _integrator1D(nullptr),
    _integrator2D(nullptr)
{
}

/**
 * @brief: Destructor.
 */
DiscreteSpaceMesh::~DiscreteSpaceMesh(void)
{
}

/**
 * @brief: Project function into space via L2.
 */
BEM::ColVector DiscreteSpaceMesh::projectFunctionL2(const ScalarFunctionBase_2D &function) const
{
    auto matrix = L2IdentityOp();
    auto rhs = testAgainstBasis(function);

    BEM::ColVector solutionV = matrix.colPivHouseholderQr().solve(rhs);
    return solutionV;
}

/**
 * @brief: Given coefficients, generate the associated function.
 * @desc: Wrapper to specific implementation on each specialization.
 */
DiscreteFunctionMeshPtr DiscreteSpaceMesh::generateFunction(const BEM::ColVector &coefficients) const
{
    std::vector<BEM::Complex> functionVector(coefficients.size(), 0);
    for (int i = 0; i < coefficients.size(); ++i) {
        functionVector[i] = coefficients[i];
    }
    return generateFunction(functionVector);
}


/**
 * @brief: Given a scalar function, test it against the basis and return the vector.
 */
BEM::ColVector DiscreteSpaceMesh::testAgainstBasis(const ScalarFunctionBase_2D &function) const
{
    auto size = getSize();
    BEM::ColVector result(size);
    for (unsigned i = 0; i < size; ++i) {
        result[i] = testAgainstBaseElement(function, i);
    }
    return result;    
}

/**
 * @brief: Constructor.
 */
DiscreteFunctionMesh::DiscreteFunctionMesh(const Mesh1D &mesh, const std::vector<BEM::Complex> coefficients) :
    _mesh{mesh},
    _coefficients{coefficients}
{
}

/**
 * @brief: returns evaluation on t. Here, t must be within the bounds of the parametrization of the discretized curve.
 */
BEM::Complex DiscreteFunctionMesh::operator()(double t) const
{
    auto location = _mesh.elementWithPoint(t);
    unsigned index = location.first;
    double position = location.second;
    return evaluate(index, position);
}

/**
 * @brief: Return L2 norm of the discrete function.
 */
double DiscreteFunctionMesh::L2Norm(void) const
{
    double result = 0.0;
    GaussLegendre_1D integrator(4); //Currently hardcoded precision. No elements higher than order 1 are implemented so 4 is more than enough.
    
    for (unsigned i = 0; i < _mesh.numElements(); ++i) {
        auto element = _mesh.getElement(i); //Iterate over the mesh elements
        ExplicitScalarFunction_1D integrand([&](double t){
            return std::abs(evaluate(i, t)*evaluate(i, t))*element.getSize();
        });
        auto localResult = integrator.integrate(0, 1, integrand); // Integrate over the element
        assert(localResult.imag() == 0.0 and "L2 local integration has imaginary part"); // Check no imaginary part.
        assert(localResult.real() >= 0.0 and "L2 local integration has negative real part"); // Check positive result
        result += localResult.real();
    }
    return std::sqrt(result);
}

/**
 * @brief: Returns L2 norm of *this - other.
 * @desc: Given a function to be evaluated over the discretized curve, this method computes the L2
 * distance between *this and the given function. The purpose is to compute errors between functions
 * and given analytical solutions to some PDE.
 */ 
double DiscreteFunctionMesh::L2Error(const ExplicitScalarFunction_2D &other) const
{
    double result = 0.0;
    GaussLegendre_1D integrator(5);
    boost::math::quadrature::tanh_sinh<double> integratorS(15);
    for (unsigned i = 0; i < _mesh.numElements(); ++i) {
        auto element = _mesh.getElement(i);
        auto integrandL = [&](double t){
            auto point = element(t);
            return std::pow(std::abs(evaluate(i, t) - other(point.getX(), point.getY())), 2.)*element.getSize();
        };
        auto localResult = integratorS.integrate(integrandL, 0, 1);
        result += std::abs(localResult);
    }
    return std::sqrt(result);
}
