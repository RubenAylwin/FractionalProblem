#ifndef GREEN_LAPLACE
#define GREEN_LAPLACE

#include <GreenLog.h>
#include <MyTypes.h>

//////////////////////////////////////
// Green function for the Laplacian //
//////////////////////////////////////

// Forward declaration
class Point2D;
class DiscreteFunction_1D;
class Curve2D;

class GreenL2D : public GreenLogSing2D {
public:
    GreenL2D(void) = default;
    BEM::Complex operator()(const Point2D &X, const Point2D &Y) const override;
private:

};

#endif
