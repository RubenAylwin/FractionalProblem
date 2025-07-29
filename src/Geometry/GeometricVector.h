#ifndef VECTOR_GEOMETRY
#define VECTOR_GEOMETRY

////////////////////////////////
// Class for geometric vector //
////////////////////////////////

class Vector2D{
public:
    Vector2D(double x, double y);
    Vector2D(const Vector2D &other) = default;
    Vector2D(Vector2D &&other) = default;

    Vector2D operator+(const Vector2D &other) const {return Vector2D(this->getX() + other.getX(), this->getY() + other.getY());}
    Vector2D operator-(const Vector2D &other) const {return Vector2D(this->getX() - other.getX(), this->getY() - other.getY());}
    Vector2D operator*(const double r) const {return Vector2D(this->getX()*r, this->getY()*r);}
    
    double norm(void) const;
    void show(void) const;

    double getX(void) const { return _x; }
    double getY(void) const { return _y; }
    
private:
    const double _x = 0.0;
    const double _y = 0.0;
};
#endif
