#include <SobolevSlobodeckii.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <DiscreteSpace.h>
#include <ScalarValuedFunction.h>
#include <Msg.h>
#include <Utilities.h>

useMessages("SOB_SLO");

/**
 * @brief: Constructor
 */
SobolevSlobodeckii::SobolevSlobodeckii(const DiscreteSpaceOnCurve_1D &space, double order):
    DifferentialOperator(space, space),
    _order(order),
    _elements()
{
    _integrator2D.reset(new GaussLegendre_2D<10, 10>());
    const auto &partition = space.getPartition();
    for (size_t p = 0; p < partition.size() - 1; ++p) {
        _elements.emplace_back(partition[p], partition[p+1]);
    }
}


/**
 * @brief: Product between to elements.
 */
BEM::Complex SobolevSlobodeckii::indexedDuality(const unsigned i, const unsigned j)
{
    auto &bsi = _trialSpace.basisFunction(i);
    auto &bsj = _trialSpace.basisFunction(j);
    auto &supportJ = bsj.brokenSupport();
    auto &supportI = bsi.brokenSupport();
    auto join = BEM::join(supportJ, supportI);
    auto Cjoin = BEM::remove(join, _elements);
    // This may not be computed with a lot of precision, since there is some singularity at t = s.
    // However, it was good enough to work since it was only needed to compute energy norms.
    // It is no longer used, since we changed to using the product between RL derivatives.
    ExplicitScalarFunction_2D integrand([&bsi, &bsj, this](double t, double s) -> BEM::Complex {
        if (t == s) {
            return 0.;
        }
        return (bsi(t) - bsi(s))*(bsj(t)-bsj(s))/std::pow((std::abs(t - s)), 1. + 2*_order);
    });
    return _integrator2D->integrate(join, join, integrand) + _integrator2D->integrate(Cjoin, join, integrand) + _integrator2D->integrate(join, Cjoin, integrand);
}
