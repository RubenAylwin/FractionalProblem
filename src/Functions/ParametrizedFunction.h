#ifndef PARAMETRIZED_FUNCTION
#define PARAMETRIZED_FUNCTION
#include <MyTypes.h>
#include <ScalarValuedFunction.h>

/////////////////////////////////////////////////////////////
// Classes for functions constructed from a parameter list //
/////////////////////////////////////////////////////////////


/**
 * @brief: Base class.
 */
class ParametrizedFunction_1D : public ScalarFunctionBase_1D {
public:
    ParametrizedFunction_1D(BEM::CVector params, double supLow, double supUpp);
    BEM::Interval1D support(void) const override;
    const BEM::Support1DL &brokenSupport(void) const override;
protected:
    BEM::CVector _params;
    BEM::Support1DL _support;
};


/**
 * @brief: Trigonometric function.
 */
class TrigonometricFunction_1D : public ParametrizedFunction_1D {
public:
    TrigonometricFunction_1D(double period, BEM::CVector params);
    BEM::Complex operator()(double t) const override;
private:
    double _period = 0.0;

};


/**
 * @brief: Polynomial.
 */
class PolynomialFunction_1D : public ParametrizedFunction_1D {
public:
    PolynomialFunction_1D(BEM::CVector params);
    BEM::Complex operator()(double t) const override;
};

/**
 * @brief: Piecewise constant.
 */
class PwConstantFunction_1D : public ParametrizedFunction_1D {
public:
    PwConstantFunction_1D(double lower, double upper, BEM::CVector params);
    BEM::Complex operator()(double t) const override;
private:
    double _lower = 0.0;
    double _upper = 0.0;
};

#endif
