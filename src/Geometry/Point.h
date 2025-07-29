#ifndef POINT_GEOMETRY
#define POINT_GEOMETRY

#include <iostream>

///////////////////////////
// Class for a 2D point. //
///////////////////////////

class Point2D{
public:
    Point2D(double x, double y);
    Point2D(const Point2D &other) = default;
    Point2D(Point2D &&other) = default;

    Point2D operator+(const Point2D &other) const {return Point2D(this->getX() + other.getX(), this->getY() + other.getY());}
    Point2D operator-(const Point2D &other) const {return Point2D(this->getX() - other.getX(), this->getY() - other.getY());}
    Point2D operator*(const double scalar) const {return Point2D(this->getX()*scalar, this->getY()*scalar);}
    bool operator==(const Point2D &other) const {return (this->_x == other.getX() and this->_y == other.getY());}
    friend std::ostream &operator<<(std::ostream &os, const Point2D &point) {return os << "(" << point._x << ", " << point._y << ")";}

    
    double norm(void) const;
    void show(void) const;

    double getX(void) const { return _x; }
    double getY(void) const { return _y; }

    void setX(double x) { _x = x; }
    void setY(double y) { _y = y; }

    static Point2D interpolate(const Point2D &first, const Point2D &second, double t);
    
private:
    double _x = 0.0;
    double _y = 0.0;
};
#endif
