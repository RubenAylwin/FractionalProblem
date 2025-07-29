#include <Mesh.h>
#include <Msg.h>
#include <Element.h>
#include <Curve.h>
#include <Point.h>
#include <algorithm>

useMessages("MESH");


/**
 * @brief: Constructor
 */
Mesh1D::Mesh1D(unsigned numElements) :
    _numElem{numElements},
    _points(numElements + 1, std::shared_ptr<Point2D>(nullptr)),
    _elements(numElements, std::shared_ptr<MeshElement1D>(nullptr))
{
}

/**
 * @brief: Get an element from its right point.
 * @input: Global number of a mesh node.
 * @output: Meshh element that has the node as its right node.
 */
const MeshElement1D &Mesh1D::elementFromBPoint(unsigned i) const
{
    unsigned index = (i >= 1) ? i - 1 : _numElem - 1;
    assert(_elements[index]->getB() == point(i));
    return *(_elements[index]);
}

/**
 * @brief: Get an element from its left point.
 * @input: Global number of a mesh node.
 * @output: Mesh element that has the node as its left node.
 */
const MeshElement1D &Mesh1D::elementFromAPoint(unsigned i) const
{
    unsigned index = (i <= _numElem) ? i : 0;
    assert(_elements[index]->getA() == point(i));
    return *(_elements[index]);
}

/**
 * @brief: Get the point in R2 associated to the ith node.
 * @input: Global number of a mesh node.
 * @output: 2D coordinates of the node.
 */
const Point2D &Mesh1D::point(unsigned i) const
{
    unsigned index = (i <= _numElem) ? i : 0;
    return *(_points[index]);
}

/**
 * @brief: Constructor.
 */
MeshCurve1D::MeshCurve1D(unsigned numElements, Curve2D &curve) :
    Mesh1D(numElements),
    _curve{curve},
    _breaks(numElements + 1, 0.0)
{
    double length = (_curve.getUppLim() - _curve.getLowLim())/numElements;
    double a = _curve.getLowLim();
    double b = a + length;
    double shift = 0.0;
    _breaks[0] = a;
    for (unsigned i = 0; i < numElements; ++i) {
        _elements[i].reset(new MeshElement1D(i, _curve.at(a), _curve.at(b), shift));
        shift += _elements[i]->getSize();
        _points[i].reset(new Point2D(_curve.at(a)));
        msg(5) << "Inserted element with points " << _curve.at(a) << " and " << _curve.at(b) << endMsg;
        a = b;
        b = a + length;
        _breaks[i+1] = a;
    }
    assert(std::abs(_breaks[numElements] - _curve.getUppLim()) < length/10.);
    _points[numElements].reset(new Point2D(_curve.at(a)));
    assert(_elements.back());
    assert(_points.back());
    assert((*(_points.back()) - _curve.at(_curve.getUppLim())).norm() < 1E-7);
    assert((*(_points.front()) - _curve.at(_curve.getLowLim())).norm() < 1E-7);
}

/**
 * @brief: Get the mesh element on which the given point lays.
 * @input: double t indicating the point on the geometry's parametrization for which we want the element.
 * @output: Pair with: <Global number of correspoding element, position of the point in the element>.
 */
std::pair<unsigned, double> MeshCurve1D::elementWithPoint(double t) const
{
    assert(t <= _curve.getUppLim() and t >= _curve.getLowLim());
    if (t == _curve.getUppLim()) {
        return std::pair<unsigned, double>(_numElem - 1u, 1.0);
    }

    if (t == _curve.getLowLim()) {
        return std::pair<unsigned, double>(0u, 0.0);
    }
    
    auto lb = std::lower_bound(_breaks.begin(), _breaks.end(), t);
    auto lowLim = *(lb - 1);
    auto uppLim = *(lb);
    assert(lowLim <= t);
    assert(uppLim >= t);
    
    return std::pair<unsigned, double>((lb - _breaks.begin()) - 1u, (t-lowLim)/(uppLim-lowLim));
}

/**
 * @brief: Constructor.
 */
MeshCurveGraded1D::MeshCurveGraded1D(unsigned numElements, Curve2D &curve, double grade) :
    Mesh1D(numElements),
    _curve{curve},
    _breaks(numElements + 1, 0.0)
{
    double totalLength = (_curve.getUppLim() - _curve.getLowLim());
    double length = totalLength/numElements;
    
    double a = _curve.getLowLim();
    double movement = std::pow(length/totalLength, grade)*totalLength;
    double b = a + movement;
    double shift = 0.0;
    _breaks[0] = a;
    for (unsigned i = 0; i < numElements; ++i) {
        _elements[i].reset(new MeshElement1D(i, _curve.at(a), _curve.at(b), shift));
        shift += _elements[i]->getSize();
        _points[i].reset(new Point2D(_curve.at(a)));
        msg(5) << "Inserted element with points " << _curve.at(a) << " and " << _curve.at(b) << endMsg;
        a = b;
        b = std::pow(length*(i+2)/totalLength, grade)*totalLength;
        _breaks[i+1] = a;
    }
    assert(std::abs(_breaks[numElements] - _curve.getUppLim()) < length/10.);
    _points[numElements].reset(new Point2D(_curve.at(a)));
    assert(_elements.back());
    assert(_points.back());
    assert((*(_points.back()) - _curve.at(_curve.getUppLim())).norm() < 1E-7);
    assert((*(_points.front()) - _curve.at(_curve.getLowLim())).norm() < 1E-7);
}

/**
 * @brief: Same as for MeshCurve1D but for this specialization.
 * TODO: Duplicated code.
 */
std::pair<unsigned, double> MeshCurveGraded1D::elementWithPoint(double t) const
{
    assert(t <= _curve.getUppLim() and t >= _curve.getLowLim());
    if (t == _curve.getUppLim()) {
        return std::pair<unsigned, double>(_numElem - 1u, 1.0);
    }

    if (t == _curve.getLowLim()) {
        return std::pair<unsigned, double>(0u, 0.0);
    }
    
    auto lb = std::lower_bound(_breaks.begin(), _breaks.end(), t);
    auto lowLim = *(lb - 1);
    auto uppLim = *(lb);
    assert(lowLim <= t);
    assert(uppLim >= t);
    
    return std::pair<unsigned, double>((lb - _breaks.begin()) - 1u, (t-lowLim)/(uppLim-lowLim));
}
