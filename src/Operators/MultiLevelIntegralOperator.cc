#include <MultiLevelIntegralOperator.h>
#include <DiscreteSpace.h>
#include <Curve.h>

/**
 * @brief: Construct the multilevel operator
 */
MultiLevelIntegralOperator::MultiLevelIntegralOperator(const DiscreteSpaceOnCurve_1D &fineTrialSpace, const DiscreteSpaceOnCurve_1D &fineTestSpace, const DiscreteSpaceOnCurve_1D &coarseTrialSpace, const DiscreteSpaceOnCurve_1D &coarseTestSpace) :
    Operator(fineTestSpace.getSize(), fineTrialSpace.getSize()),
    _fineTrialSpace(fineTrialSpace),
    _fineTestSpace(fineTestSpace),
    _coarseTrialSpace(coarseTrialSpace),
    _coarseTestSpace(coarseTestSpace)
{
}

/**
 * @brief: Explicit destructor
 */
MultiLevelIntegralOperator::~MultiLevelIntegralOperator()
{
}
