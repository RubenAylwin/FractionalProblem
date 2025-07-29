#ifndef GREEN_FUNCTION_BASE
#define GREEN_FUNCTION_BASE

#include <MyTypes.h>
#include <vector>
#include <cassert>

/////////////////////////////////////////////////
// Base class for green functions on 2D space. //
/////////////////////////////////////////////////

//Forward declarations
class BasisFunction_1D;
class Curve2D;
class Point2D;
class MeshElement1D;
class BasisFunctionMesh;

class GreenFunction2D {
public:
    GreenFunction2D(void) {};
    virtual ~GreenFunction2D(void) = default;
    virtual BEM::Complex operator()(const Point2D &X, const Point2D &Y) const = 0;
    virtual BEM::Complex singularity(const Point2D &X, const Point2D &Y) const = 0;
    virtual BEM::Complex singularitySimple(const Point2D &X, const Point2D &Y) const = 0;
    virtual BEM::Complex integrateSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction) const = 0;
    virtual BEM::Complex integrateSingularity([[maybe_unused]] const MeshElement1D &testElement, [[maybe_unused]] const MeshElement1D &trialElement, [[maybe_unused]] const BasisFunctionMesh &testFunction, [[maybe_unused]] const BasisFunctionMesh &trialFunction) const {return 0.0;}
    virtual BEM::Complex smoothPart(const Point2D &X, const Point2D &Y) const;
private:
};

#endif
