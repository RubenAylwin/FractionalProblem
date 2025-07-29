#include <ParametrizedFunction.h>
#include <Msg.h>
useMessages("PAR_FUN");

/**
 * @brief: Constructor. Just saves the parameters.
 */
ParametrizedFunction_1D::ParametrizedFunction_1D(BEM::CVector params, double supLow, double supUpp) :
    _params{params},
    _support{BEM::Interval1D(supLow, supUpp)}
{
    msg(5) << "start::ParametrizedFunction_1D::ParametrizedFunction_1D()" << endMsg;
    msg(5) << "end::ParametrizedFunction_1D::ParametrizedFunction_1D()" << endMsg;
}

/**
 * @brief: Return the support. Needed for completion.
 */
BEM::Interval1D ParametrizedFunction_1D::support(void) const
{
    return _support.front();
}

/**
 * @brief: Return the support. Needed for completion.
 */
const BEM::Support1DL &ParametrizedFunction_1D::brokenSupport(void) const
{
    return _support;
}

/**
 * @brief: Constructor. Saves parameters and period.
 */
TrigonometricFunction_1D::TrigonometricFunction_1D(double period, BEM::CVector params) :
    ParametrizedFunction_1D(params, 0., period),
    _period{period}
{
    msg(5) << "start::TrigonometricFunction_1D::TrigonometricFunction_1D()" << endMsg;
    msg(5) << "end::TrigonometricFunction_1D::TrigonometricFunction_1D()" << endMsg;
}

/**
 * @brief: Evaluation.
 */
BEM::Complex TrigonometricFunction_1D::operator()(double t) const
{
    BEM::Complex result = 0.0;
    int counter = .0;//To keep track of the order.
    for (size_t i = 0u; i < _params.size(); ++i) {
        auto param = _params[i];
        //Odd parameters correspond to Cos. Even correspond to sine.
        result += (i%2 == 0) ? param*std::cos(2.*M_PI*t*counter/_period) : param*std::sin(2.*M_PI*t*counter/_period);
        if (i%2==0) {
            ++counter;
        }
    }
    return result;
}

/**
 * @brief: Constructor. Saves parameters.
 */
PolynomialFunction_1D::PolynomialFunction_1D(BEM::CVector params) :
    ParametrizedFunction_1D(params, -1000, 1000)
{
    msg(5) << "start::PolynomialFunction_1D::PolynomialFunction_1D()" << endMsg;
    msg(5) << "end::PolynomialFunction_1D::PolynomialFunction_1D()" << endMsg;
}

/**
 * @brief: Evaluation.
 */
BEM::Complex PolynomialFunction_1D::operator()(double t) const
{
    BEM::Complex result = 0.0;
    for (size_t i = 0u; i < _params.size(); ++i) {
        result += _params[i]*std::pow(t, i);
    }
    return result;
}

/**
 * @brief: Constructor. Saves parameters.
 */
PwConstantFunction_1D::PwConstantFunction_1D(double lower, double upper, BEM::CVector params) :
    ParametrizedFunction_1D(params, lower, upper),
    _lower{lower},
    _upper{upper}
{
    msg(5) << "start::PwConstantFunction_1D::PwConstantFunction_1D()" << endMsg;
    msg(5) << "end::PwConstantFunction_1D::PwConstantFunction_1D()" << endMsg;
}

/**
 * @brief: Evaluation.
 */
BEM::Complex PwConstantFunction_1D::operator()(double t) const
{
    if (t < _lower or t > _upper) {
        return 0.0;
    }
    double delta = (_upper - _lower)/_params.size();
    unsigned index = std::floor((t - _lower)/delta);
    
    return _params[index];
}
