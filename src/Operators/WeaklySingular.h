#ifndef WEAKLY_SINGULAR_OPERATOR
#define WEAKLY_SINGULAR_OPERATOR

#include <MyTypes.h>
#include <IntegralOperator.h>
#include <vector>
#include <mutex>
#include <memory>


////////////////////////////////////////////////
// Classes for Weakly Singular Operator (BEM) //
////////////////////////////////////////////////

// Forward declarations
class GreenFunction2D;
class DiscreteSpaceOnCurve_1D;
class DiscreteSpaceMesh;
class ScalarFunctionBase_2D;
class ExplicitScalarFunction_2D;
class Transformation_2D;

/**
 * @brief: Weakly singular operator on Discrete spaces.
 * DEPRECATED in favour of mesh version.
 */
class WeaklySingular : public IntegralOperator {
public:
    WeaklySingular(const DiscreteSpaceOnCurve_1D &trialSpace, const DiscreteSpaceOnCurve_1D &testSpace, const GreenFunction2D &greenFunc);
    ~WeaklySingular(void);
    void setHighPrecision(void);
    const GreenFunction2D &getGreen();
    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
    
private:
    const GreenFunction2D &_greenFunc;
    
};

/**
 * @brief: Weakly singular operator on discrete spaces defined over a mesh.
 */
class WeaklySingularMesh : public IntegralOperatorMesh {
public:
    WeaklySingularMesh(const DiscreteSpaceMesh &trialSpace, const DiscreteSpaceMesh &testSpace, const GreenFunction2D &greenFunc);
    WeaklySingularMesh(const DiscreteSpaceMesh &trialSpace, const DiscreteSpaceMesh &testSpace, const GreenFunction2D &greenFunc, Transformation_2D &&transformation, ExplicitScalarFunction_2D &&jacobian);
    ~WeaklySingularMesh(void);
    void setHighPrecision(void);
    const GreenFunction2D &getGreen();
    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
    
private:
    const GreenFunction2D &_greenFunc;
    std::unique_ptr<Transformation_2D> _transformation = nullptr;
    std::unique_ptr<ExplicitScalarFunction_2D> _jacobian = nullptr;    
};

#endif
