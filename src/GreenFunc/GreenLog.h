#ifndef GREEN_LOG_SING
#define GREEN_LOG_SING

#include <Green.h>
#include <MyTypes.h>
#include <GreenQP.h>

//////////////////////////////////////////////////////////////
// Base class for functions with a logarithmic singularity. //
//////////////////////////////////////////////////////////////

//Forward declarations
class BasisFunction_1D;
class Point2D;
class Curve2D;
class Integrator_1D;
class Integrator_2D;

class GreenLogSing2D : public PeriodizableGreenFunction2D {
public:
    GreenLogSing2D(void);
    ~GreenLogSing2D(void);
    virtual void set1DQuad(unsigned N);
    virtual void setHighPrecision(int level = 0);
    virtual BEM::Complex singularity(const Point2D &X, const Point2D &Y) const override;
    virtual BEM::Complex singularitySimple(const Point2D &X, const Point2D &Y) const override;
    BEM::Complex integrateSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction) const override;
    BEM::Complex integrateSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction) const override;
    BEM::Complex integrateShiftedSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction, const double shift) const override;
    BEM::Complex integrateShiftedSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction, const double shift) const override;
protected:
    mutable int quadPoints = 2;
    mutable std::unique_ptr<Integrator_1D> _integrator1D;
    mutable std::unique_ptr<Integrator_2D> _integrator2D;
};

#endif
