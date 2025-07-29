#ifndef ELEMENT
#define ELEMENT
#include <MyTypes.h>
#include <Point.h>
#include <string>

/**
 * @brief: Element of a partition. DEPRECATED in favour of MeshElement1D and meshes instead of partitions.
 */
class Element_1D {
 public:
    Element_1D(unsigned a, unsigned b, unsigned number, unsigned dof) : _a{a}, _b{b}, _number{number}, _dof{dof} {};
    unsigned getA(void) const {return _a;}
    unsigned getB(void) const {return _b;}
    unsigned getNum(void) const {return _number;}
    unsigned getDof(void) const {return _dof;}
    bool operator==(const Element_1D &other) const {return (this->_a == other.getA()) and (this->_b == other.getB());}
 private:
    unsigned _a;
    unsigned _b;
    unsigned _number;
    unsigned _dof;
};

/**
 * @brief: Element of a mesh.
 */
class MeshElement1D {
 public:
    MeshElement1D(unsigned index, Point2D a, Point2D b, double shift) : _index{index}, _a{a}, _b{b}, _shift{shift} {}
;
    Point2D getA(void) const {return _a;}
    Point2D getB(void) const {return _b;}
    double getSize(void) const {return (_b - _a).norm();}
    unsigned getIndex(void) const {return _index;}
    Point2D operator()(double t) const {return _a + (_b-_a)*t;}
    double operator[](double t) const {return _shift+t*getSize();}
    bool operator==(const MeshElement1D &other) const {return (this->_a == other.getA()) and (this->_b == other.getB());}
    friend std::ostream &operator<<(std::ostream &os, const MeshElement1D &element) {return os << "[" << element._a << ", " << element._b << "]";}
 private:
    unsigned _index;
    Point2D _a;
    Point2D _b;
    double _shift = 0.0;
};


#endif
