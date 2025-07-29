#ifndef CURVE_GEOMETRY
#define CURVE_GEOMETRY

#include <vector>
#include <Point.h>
#include <GeometricVector.h>
#include <memory>

///////////////////////////////////////////////////////////////////////////////
// Classes representing different types of 1-D curves embedded on 2-D space. //
///////////////////////////////////////////////////////////////////////////////

//Forward declaration
class ScalarFunctionBase_1D;

/**
 * @brief: Base class.
 */
class Curve2D {
public:
    Curve2D(double lowLimit, double uppLimit) : _lowLimit{lowLimit}, _uppLimit{uppLimit} {}
    virtual ~Curve2D(void) = default;
    virtual Point2D at(const double t) const;
    virtual Vector2D normal(const double t) const;
    virtual double jacobian(const double t) const;
    virtual std::unique_ptr<ScalarFunctionBase_1D> surfaceMeasure(void) const;
    virtual Point2D evaluateAt(const double t) const = 0;
    virtual Vector2D evaluateNormalAt(const double t) const = 0;
    virtual std::vector<double> getParameters(void) const {return std::vector<double>{};}
    virtual double getLowLim(void) {return _lowLimit;}
    virtual double getUppLim(void) {return _uppLimit;}
private:
    double _lowLimit;
    double _uppLimit;
};

/**
 * @brief: Curve describing a straight segment.
 */
class StraightCurve : public Curve2D {
public:
    StraightCurve(double length);
    Point2D evaluateAt(const double t) const override;
    Vector2D evaluateNormalAt(const double t) const override;
private:
    double _length;
};

/**
 * @brief: A periodic curve.
 */
class PeriodicCurve : public Curve2D {
public:
    PeriodicCurve(double period) : Curve2D(0, 1), _period{period} {}
    double getPeriod(void) {return _period;}
protected:
    double _period = 0;
};

/**
 * @brief: Specialization of periodic curve to curves given by trigonometric functions.
 */
class TrigonometricCurve : public PeriodicCurve {
public:
    TrigonometricCurve(double period, double height, std::vector<double> sineCoefficents, std::vector<double> cosineCoefficents);
    std::vector<double> getParameters(void) const override;
private:
    Point2D evaluateAt(const double t) const override;
    Vector2D evaluateNormalAt(const double t) const override;
    const double _height;
    const std::vector<double> _sineCoefficients;
    const std::vector<double> _cosineCoefficients;
};

/**
 * @brief: Specialization of periodic curve to a periodic polygonal curve.
 */
class PolyPeriodicCurve : public PeriodicCurve {
public:
    PolyPeriodicCurve(double period, std::vector<Point2D> &&orderedPoints);
    std::vector<double> getPartition(void) {return _intervalPartition;}
private:
    Point2D evaluateAt(const double t) const override;
    Vector2D evaluateNormalAt(const double t) const override;
    const unsigned _numLines;
    std::vector<Point2D> _orderedPoints;
    std::vector<double> _intervalPartition;
};

#endif
