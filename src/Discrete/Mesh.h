#ifndef MESH
#define MESH

#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cassert>

/////////////////////////////////////////////////////////////////////////////////////////
// This file contains all clases for meshes on 1D curves (which may be embedded on 2D) //
/////////////////////////////////////////////////////////////////////////////////////////

// Forward Declarations
class MeshElement1D;
class Curve2D;
class Point2D;

/**
 * @brief: Mesh class (base class).
 */
class Mesh1D {
public:
    virtual ~Mesh1D(void) {};
    virtual const MeshElement1D &elementFromBPoint(unsigned i) const;
    virtual const MeshElement1D &elementFromAPoint(unsigned i) const;
    virtual const MeshElement1D &getElement(unsigned i) const {assert(i < _numElem); return *_elements[i];}
    virtual const Point2D &point(unsigned i) const;
    virtual std::pair<unsigned, double> elementWithPoint(double t) const = 0;
    virtual unsigned numElements(void) const {return _numElem;};
protected:
    Mesh1D(unsigned numElements);
    unsigned _numElem = 0U;
    std::vector<std::shared_ptr<Point2D>> _points;
    std::vector<std::shared_ptr<MeshElement1D>> _elements;
};


/**
 * @brief: Uniform mesh over a curve. Uniformity is over given parametrization.
 */
class MeshCurve1D : public Mesh1D {
public:
    MeshCurve1D(unsigned numElements, Curve2D &curve);
    std::pair<unsigned, double> elementWithPoint(double t) const override;
private:
    Curve2D &_curve;
    std::vector<double> _breaks;
};


/**
 * @brief: Graded mesh over a curve. Gradation is over given parametrization.
 */
class MeshCurveGraded1D : public Mesh1D {
public:
    MeshCurveGraded1D(unsigned numElements, Curve2D &curve, double grade);
    std::pair<unsigned, double> elementWithPoint(double t) const override;
private:
    Curve2D &_curve;
    std::vector<double> _breaks;
};

#endif
