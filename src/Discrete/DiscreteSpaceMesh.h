#ifndef DISCRETE_MESH_SPACE
#define DISCRETE_MESH_SPACE

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

////////////////////////////////////////////////////////////////////////////////////////////////////
// This file contains declarations for the base clases of discrete spaces defined over 1D meshes. //
////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Mesh1D;
class Integrator_1D;
class Integrator_2D;

//Discrete Function of R.
class DiscreteFunctionMesh;
using DiscreteFunctionMeshPtr = std::unique_ptr<DiscreteFunctionMesh>;

/**
 * @brief: Discrete function over a mesh (base class).
 */
class DiscreteFunctionMesh {
public:
    DiscreteFunctionMesh(const Mesh1D &mesh, const std::vector<BEM::Complex> coefficients);
    DiscreteFunctionMesh(void) = delete;
    DiscreteFunctionMesh(const DiscreteFunctionMesh &other) = delete;
    virtual const Mesh1D &getMesh(void) const {return _mesh;}
    virtual const std::vector<unsigned> &support(void) const {return _support;}
    virtual ~DiscreteFunctionMesh(void) {};
    virtual BEM::Complex evaluate(unsigned elementIndex, double t) const = 0;
    virtual BEM::Complex operator()(double t) const;
    virtual std::vector<BEM::Complex> getCoefficients(void) {return _coefficients;}
    virtual double L2Norm(void) const;
    virtual double L2Error(const ExplicitScalarFunction_2D &other) const;
    virtual void operator-=(DiscreteFunctionMesh &other) {_coefficients = _coefficients - other._coefficients;};
    virtual void operator*=(DiscreteFunctionMesh &other) {
        for (size_t i = 0; i < _coefficients.size(); ++i) {
            _coefficients[i] = _coefficients[i] * other._coefficients[i];
        }};
    virtual DiscreteFunctionMesh &derivative(void) const {
        assert(false and "Not implemented yet.");
    };
    virtual DiscreteFunctionMesh &derivative([[maybe_unused]] int order) const {
        assert(false and "Not implemented yet.");
        return derivative();};
protected:
    const Mesh1D &_mesh;
    std::vector<unsigned> _support;
    std::vector<BEM::Complex> _coefficients;

};

/**
 * @brief: Basis function on a Mesh (base class).
 */
class BasisFunctionMesh : virtual public DiscreteFunctionMesh
{
public:
    BasisFunctionMesh(const Mesh1D &mesh, const std::vector<BEM::Complex> coefficients) : DiscreteFunctionMesh(mesh, coefficients) {};
    BasisFunctionMesh(void) = delete;
    virtual ~BasisFunctionMesh(void) {};
    virtual const ScalarFunctionBase_1D &antiDerivative(unsigned elementIndex) const = 0;
protected:
    mutable std::unordered_map<unsigned, ScalarFunctionBasePtr_1D> _antiDerivatives;
};

// Typenames
using DiscreteFunctionMeshPtr = std::unique_ptr<DiscreteFunctionMesh>;
using DiscreteFunctionalMesh = std::function<BEM::Complex (const ScalarFunctionBase_1D &function)>;

/**
 * @brief: Discrete space over a mesh (base class).
 */
class DiscreteSpaceMesh {
public:
    DiscreteSpaceMesh(const Mesh1D &mesh);
    virtual ~DiscreteSpaceMesh(void);
    // Test Function against Space
    virtual BEM::ColVector testAgainstBasis(const ScalarFunctionBase_2D &function) const;
    virtual BEM::Complex testAgainstBaseElement(const ScalarFunctionBase_2D &function, unsigned globalNumber) const = 0;
    // Get Discrete Functions
    virtual DiscreteFunctionMeshPtr generateFunction(const BEM::ColVector &coefficients) const;
    virtual DiscreteFunctionMeshPtr generateFunction(const std::vector<BEM::Complex> coefficients) const = 0;
    virtual DiscreteFunctionMeshPtr generateFunction([[maybe_unused]] const ScalarFunctionBase_2D &function) const {
        assert(false and "not implemented");
        return nullptr;};
    virtual BEM::ColVector projectFunctionL2(const ScalarFunctionBase_2D &function) const;
    // FEM specifics
    virtual const BasisFunctionMesh &basisFunction(const unsigned globalNumber) const = 0;
    virtual const Mesh1D &getMesh(void) const {return _mesh;}
    // Info
    virtual unsigned getSize(void) const = 0;
    virtual const std::string &getBaseName(void) const = 0;
    // Aux
    virtual BEM::Matrix L2IdentityOp(void) const = 0;
    
protected:
    const Mesh1D &_mesh;
    std::shared_ptr<Integrator_1D> _integrator1D;
    std::shared_ptr<Integrator_2D> _integrator2D;
};

#endif
