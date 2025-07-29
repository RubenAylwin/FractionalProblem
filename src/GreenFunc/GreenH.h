#ifndef GREEN_HELMHOLTZ
#define GREEN_HELMHOLTZ

#include <GreenLog.h>
#include <MyTypes.h>

/////////////////////////////////////////
// Class for Helmholtz green function. //
/////////////////////////////////////////

// Forward declarations
class Point2D;
class DiscreteFunction_1D;
class Curve2D;

class GreenH2D : public GreenLogSing2D {
public:
    GreenH2D(double wavenumber);
    BEM::Complex operator()(const Point2D &X, const Point2D &Y) const override;
    double getWavenumber(void) { return _wavenumber; };
private:
    const double _wavenumber = 1.0;
};

#endif
    
