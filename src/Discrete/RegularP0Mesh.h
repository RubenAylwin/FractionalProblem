#ifndef REGULAR_P0_MESH
#define REGULAR_P0_MESH

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
#include <DiscreteSpaceMesh.h>

//////////////////////////////////////////////////////////////////////////
// Class for a regular space of piecewise constant functions on a mesh. //
//////////////////////////////////////////////////////////////////////////

// Forward declarations
class Integrator_1D;
class Integrator_2D;

/**
 * @brief: Class for piecewise P1 functions on a mesh.
 */
class RegularP0Mesh_1D : public DiscreteSpaceMesh {
public:
    RegularP0Mesh_1D(const Mesh1D &mesh);
    ~RegularP0Mesh_1D(void);
    // Test Function against Space
    BEM::Complex testAgainstBaseElement(const ScalarFunctionBase_2D &function, unsigned globalNumber) const override;
    // Get Discrete Function
    DiscreteFunctionMeshPtr generateFunction(const std::vector<BEM::Complex> coefficients) const override;
    DiscreteFunctionMeshPtr generateFunction(const ScalarFunctionBase_2D &function) const override;
    // FEM Specifics
    const BasisFunctionMesh &basisFunction(const unsigned globalNumber) const override;
    // Info
    unsigned getSize(void) const override { return _size; };
    const std::string &getBaseName(void) const override { return _baseName; };

    // Aux
    virtual BEM::Matrix L2IdentityOp(void) const override { return BEM::Matrix(1,1);}

    // Basis
    class P0Function : virtual public DiscreteFunctionMesh {
    public:
        P0Function(void) = delete;
        BEM::Complex evaluate(unsigned elementIndex, double t) const override;
    private:
        P0Function(const Mesh1D &_mesh, const std::vector<BEM::Complex> coefficients);
        DiscreteFunctionMesh &derivative(void) const override;
        friend class RegularP0Mesh_1D;
        friend class RegularP1_0Mesh_1D;
        friend class RegularP1_Mesh_1D;
        mutable std::mutex _derMutex;
        mutable std::mutex _supportMutex;
        mutable std::unique_ptr<DiscreteFunctionMesh> _derivative = nullptr;
    };

    class BasisFunction : public BasisFunctionMesh, public P0Function {
    public:
        // BEM::Complex operator()(double t) const override;
    private:
        BasisFunction(const Mesh1D &mesh, unsigned index);
        BasisFunction(void) = delete;
        const ScalarFunctionBase_1D &antiDerivative(unsigned elementIndex) const override;
        mutable std::unique_ptr<ScalarFunctionBase_2D> _antiDerivative = nullptr;

        friend class RegularP0Mesh_1D;
    };

private:
    const unsigned _size;
    const std::string _baseName = "RegularP0Mesh";
    mutable std::unordered_map<int, std::unique_ptr<BasisFunctionMesh>> _basis;
    mutable std::unordered_map<int, std::unique_ptr<DiscreteFunctionalMesh>> _dof;
    mutable std::mutex _basisMutex;
    mutable std::mutex _dofMutex;
};

#endif
