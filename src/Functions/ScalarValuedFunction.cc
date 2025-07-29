#include <ScalarValuedFunction.h>
#include <Curve.h>
#include <iostream>
#include <list>

/**
 * @brief: Product between two functions.
 */
ScalarFunctionBasePtr_1D ScalarFunctionBase_1D::operator*(const ScalarFunctionBase_1D &other) const
{
    auto product = [&](const double t) { return (*this)(t)*other(t); };
    return ScalarFunctionBasePtr_1D(new ExplicitScalarFunction_1D(std::move(product)));
}

/**
 * @brief: Product between function and scalar.
 */
ScalarFunctionBasePtr_1D ScalarFunctionBase_1D::operator*(const BEM::Complex scalar) const
{
    auto product = [&](const double t) { return (*this)(t)*scalar; };
    return ScalarFunctionBasePtr_1D(new ExplicitScalarFunction_1D(std::move(product)));
}

/**
 * @brief: Addition between two functions.
 */
ScalarFunctionBasePtr_1D ScalarFunctionBase_1D::operator+(const ScalarFunctionBase_1D &other) const
{
    auto sum = [&](const double t) { return (*this)(t)+other(t); };
    return ScalarFunctionBasePtr_1D(new ExplicitScalarFunction_1D(std::move(sum)));
}

/**
 * @brief: Destructor.
 */
ScalarFunctionBase_1D::~ScalarFunctionBase_1D(void)
{
}

/**
 * @brief: Tensorize two 1-D functions and get a 2-D function.
 */
ScalarFunctionBasePtr_2D ScalarFunctionBase_1D::tensor(const ScalarFunctionBase_1D &first, const ScalarFunctionBase_1D &second)
{
    auto product = [&](const double t, const double s) {return first(t)*second(s);};
    return ScalarFunctionBasePtr_2D(new ExplicitScalarFunction_2D(std::move(product)));
}


/**
 * @brief: Evaluate on t.
 */
BEM::Complex ExplicitScalarFunction_1D::operator()(double t) const
{
    return _function(t);
}

/**
 * @brief: Get the support.
 */
BEM::Interval1D ExplicitScalarFunction_1D::support(void) const
{
    return BEM::Interval1D(std::numeric_limits<double>::min(), std::numeric_limits<double>::max());
}

/**
 * @brief: Constructor.
 */
BoundaryScalarTrace_1D::BoundaryScalarTrace_1D(BEM::OneDimFunction &&function, const Curve2D &curve) :
    ExplicitScalarFunction_1D(std::forward<BEM::OneDimFunction>(function)),
     _curve{curve}
{
}

/**
 * @brief: Constructor.
 */
BoundaryScalarRestriction_1D::BoundaryScalarRestriction_1D(BEM::TwoDimFunction &&function, const Curve2D &curve) :
    BoundaryScalarTrace_1D([&](const double t) {return _function(curve.at(t).getX(), curve.at(t).getY());}, curve),
    _function(function)
{
}

/**
 * @brief: Product between 2-D functions.
 */
ScalarFunctionBasePtr_2D ScalarFunctionBase_2D::operator*(const ScalarFunctionBase_2D &other) const
{
    auto product = [&](const double t, const double s) { return (*this)(t, s)*other(t, s); };
    return ScalarFunctionBasePtr_2D(new ExplicitScalarFunction_2D(std::move(product)));
}

/**
 * @brief: Evaluation.
 */
BEM::Complex ExplicitScalarFunction_2D::operator()(const double t, const double s) const
{
    return _function(t, s);
}
