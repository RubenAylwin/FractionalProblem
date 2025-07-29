#ifndef GREEN_QUASI_PERIODIC
#define GREEN_QUASI_PERIODIC

#include <Green.h>
#include <MyTypes.h>
#include <memory>

////////////////////////////////////////////////////////////////////////////////
// Base class and other classes for periodizable and periodic Green functions //
////////////////////////////////////////////////////////////////////////////////

/**
 * @brief: Base class for green functions that are periodic
 */
class PeriodizableGreenFunction2D : public GreenFunction2D {
public:
    PeriodizableGreenFunction2D(void) {};
    virtual BEM::Complex integrateShiftedSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction, const double shift) const = 0;
    virtual BEM::Complex integrateShiftedSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction, const double shift) const = 0;
};

/**
 * @brief: Base class for quasi-periodic green functions
 */
class GreenQP2D: public GreenFunction2D {
public:
    GreenQP2D(double period);
    ~GreenQP2D(void);
    BEM::Complex operator()(const Point2D &X, const Point2D &Y) const override;
    virtual BEM::Complex singularity(const Point2D &X, const Point2D &Y) const override;
    virtual BEM::Complex singularitySimple(const Point2D &X, const Point2D &Y) const override;
    virtual BEM::Complex displacedSum(const Point2D &X, const Point2D &Y) const;
    virtual BEM::Complex windowedSum(const Point2D &X, const Point2D &Y) const;
    BEM::Complex integrateSingularity(const Curve2D &curve, const BasisFunction_1D &testFunction, const BasisFunction_1D &trialFunction) const override;
    BEM::Complex integrateSingularity(const MeshElement1D &testElement, const MeshElement1D &trialElement, const BasisFunctionMesh &testFunction, const BasisFunctionMesh &trialFunction) const override;
    virtual BEM::Complex getQP(void) { return _qp; };
    virtual void setWindowTerms(int terms) { _windowTerms = terms; };
protected:
    BEM::Complex _qp = 1.0;
    std::unique_ptr<PeriodizableGreenFunction2D> _green;
    double _period = 1.0;
    int _windowTerms = 20;
};

/**
 * @brief: Laplace QP green function
 */
class GreenLQP2D : public GreenQP2D {
public:
    GreenLQP2D(double period);
private:
};

/**
 * @brief: Helmholtz QP green function
 */
class GreenH2D;
class GreenHQP2D : public GreenQP2D {
public:
    GreenHQP2D(double period, double angle, double wavenumber);
    GreenHQP2D(double period, double angle, double wavenumber, unsigned quadPoints1D);
    BEM::Complex spectralSum(const Point2D &X, const Point2D &Y) const;
    void setHighPrecision(int i=0);
private:
    const double _angle = 0.0;
    const double _wavenumber = 1.0;
    GreenH2D *_greenH;
};

#endif

