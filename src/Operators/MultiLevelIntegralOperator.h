#ifndef MULTI_LEVEL_OPERATOR
#define MULTI_LEVEL_OPERATOR

#include <IntegralOperator.h>
#include <Identity.h>
#include <MyTypes.h>
#include <Utilities.h>
#include <DiscreteSpace.h>
#include <DiscreteSpaceMatrixMgr.h>
#include <memory>
#include <iostream>
class Integrator_1D;
class Integrator_2D;
class GreenFunction2D;

////////////////////////////////////////////////////////////////////
// Classes for multilevel operators (currently under development) //
////////////////////////////////////////////////////////////////////

/**
 * @brief: Base class for multilevel operator.
 */
class MultiLevelIntegralOperator : public Operator
{
public:
    MultiLevelIntegralOperator(const DiscreteSpaceOnCurve_1D &fineTrialSpace, const DiscreteSpaceOnCurve_1D &fineTestSpace, const DiscreteSpaceOnCurve_1D &coarseTrialSpace, const DiscreteSpaceOnCurve_1D &coarseTestSpace);
    virtual ~MultiLevelIntegralOperator();
    
protected:
    const DiscreteSpaceOnCurve_1D &_fineTrialSpace;
    const DiscreteSpaceOnCurve_1D &_fineTestSpace;
    const DiscreteSpaceOnCurve_1D &_coarseTrialSpace;
    const DiscreteSpaceOnCurve_1D &_coarseTestSpace;
};

/**
 * @brief: Specific for BEM operators. Templetized on the specific operator.
 */
template <class OperatorType>
class MultiLevelBEMOperator : public MultiLevelIntegralOperator
{
public:
    MultiLevelBEMOperator(const DiscreteSpaceOnCurve_1D &fineTrialSpace, const DiscreteSpaceOnCurve_1D &fineTestSpace, const DiscreteSpaceOnCurve_1D &coarseTrialSpace, const DiscreteSpaceOnCurve_1D &coarseTestSpace, const GreenFunction2D &greenFunc);
    BEM::Complex indexedDuality(const unsigned i, const unsigned j) override;
    void assembleMassMatrix(void) override;
    const BEM::Matrix &getFineMatrix(void) { return _fineOperator.getMatrix(); }
    const BEM::Matrix &getCoarseMatrix(void) { return _coarseOperator.getMatrix(); }
    double error(void);
private:
    OperatorType _fineOperator;
    OperatorType _coarseOperator;
    Identity _coarseTrial;
    Identity _coarseTest;
    Identity _fineToCoarseTrial;
    Identity _fineToCoarseTest;
};

/**
 * @brief: Constructor.
 */
template <class OperatorType>
MultiLevelBEMOperator<OperatorType>::MultiLevelBEMOperator(const DiscreteSpaceOnCurve_1D &fineTrialSpace, const DiscreteSpaceOnCurve_1D &fineTestSpace, const DiscreteSpaceOnCurve_1D &coarseTrialSpace, const DiscreteSpaceOnCurve_1D &coarseTestSpace, const GreenFunction2D &greenFunc) :
    MultiLevelIntegralOperator(fineTrialSpace, fineTestSpace, coarseTrialSpace, coarseTestSpace),
    _fineOperator(_fineTrialSpace, _fineTestSpace, greenFunc),
    _coarseOperator(_coarseTrialSpace, _coarseTestSpace, greenFunc),
    _coarseTrial(_coarseTrialSpace, _coarseTrialSpace),
    _coarseTest(_coarseTestSpace, _coarseTestSpace),
    _fineToCoarseTrial(_fineTrialSpace, _coarseTrialSpace),
    _fineToCoarseTest(_fineTestSpace, _coarseTestSpace)
{
}

/**
 * @brief: Assembling.
 */
template <class OperatorType>
void MultiLevelBEMOperator<OperatorType>::assembleMassMatrix(void)
{
    _fineOperator.assembleMassMatrix();

    // Matrix declarations
    //BEM::Matrix LHCF(_fineTestSpace.getSize(), _coarseTestSpace.getSize());
    //BEM::Matrix RHCF(_coarseTrialSpace.getSize(), _fineTrialSpace.getSize());
    //BEM::Matrix LHFC(_coarseTestSpace.getSize(), _fineTestSpace.getSize());
    //BEM::Matrix RHFC(_fineTrialSpace.getSize(), _coarseTrialSpace.getSize());
    
    ////////////////////////////////////////////////////////////////////////////////////
    // This (discarded) block of code computes the DoF matrices from the given spaces //
    ////////////////////////////////////////////////////////////////////////////////////
#if 0
    for (int row = 0; row < LHCF.rows(); ++row) {
        auto basis = _fineTestSpace.basisFunction(row);
        for (int col = 0; col < LHCF.cols(); ++col) {
            auto DoF = _coarseTestSpace.degreeOfFreedom(col);
            LHCF(row, col) = DoF(*basis);
        }
    }

    for (int col = 0; col < RHCF.cols(); ++col) {
        auto basis = _fineTrialSpace.basisFunction(col);
        for (int row = 0; row < RHCF.rows(); ++row) {
            auto DoF = _coarseTrialSpace.degreeOfFreedom(row);
            RHCF(row, col) = DoF(*basis);
        }
    }
    
    for (int row = 0; row < LHFC.rows(); ++row) {
        auto basis = _coarseTestSpace.basisFunction(row);
        for (int col = 0; col < LHFC.cols(); ++col) {
            auto DoF = _fineTestSpace.degreeOfFreedom(col);
            LHFC(row, col) = DoF(*basis);
        }
    }

    for (int col = 0; col < RHFC.cols(); ++col) {
        auto basis = _coarseTrialSpace.basisFunction(col);
        for (int row = 0; row < RHFC.rows(); ++row) {
            auto DoF = _fineTrialSpace.degreeOfFreedom(row);
            RHFC(row, col) = DoF(*basis);
        }
    }
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////
    // The matrices above depend only on the polinomial space and not on the geometry, so we //
    // implement a singleton that saves the matrices to avoid unnecessary computations       //
    ///////////////////////////////////////////////////////////////////////////////////////////
    auto &matrixMgr = DiscreteSpaceMatrixMgr_1D::get();
    auto &LHCF = matrixMgr.getData(_coarseTestSpace, _fineTestSpace).transpose();
    auto &RHCF = matrixMgr.getData(_coarseTrialSpace, _fineTrialSpace);    
    ////////////////////////////////////////////////////////////////////////////////////////
    // This block of code computes the multilevel matrix by computing the matrix from the //
    // coarse level and extending it to the fine level.                                   //
    ////////////////////////////////////////////////////////////////////////////////////////
#if 1
    _coarseOperator.assembleMassMatrix();
    _massMatrix.reset(new BEM::Matrix(_fineOperator.getMatrix() - BEM::product(BEM::product(LHCF, _coarseOperator.getMatrix()), RHCF)));
    return;
#endif
    ///////////////////////////////////////////////////////////////////////////////////////
    // The remaining code computes the difference operator directly from the fine matrix //
    ///////////////////////////////////////////////////////////////////////////////////////
    auto &LHFC = matrixMgr.getData(_fineTestSpace, _coarseTestSpace).transpose();
    auto &RHFC = matrixMgr.getData(_fineTrialSpace, _coarseTrialSpace);
    auto LH = BEM::product(LHCF, LHFC);
    auto RH = BEM::product(RHFC, RHCF);
    // std::cout << RH << std::endl;
    _massMatrix.reset(new BEM::Matrix(_fineOperator.getMatrix() - BEM::product(BEM::product(LH, _fineOperator.getMatrix()), RH)));
    
}


/**
 * @brief: Efficient computation of singular matrix terms for later approximation via MEIM.
 */
template <class OperatorType>
BEM::Complex MultiLevelBEMOperator<OperatorType>::indexedDuality(const unsigned i, const unsigned j)
{
    auto fineData = _fineOperator.indexedDuality(i,j);
    BEM::Matrix LHCF(_fineTestSpace.getSize(), _coarseTestSpace.getSize());
    BEM::Matrix RHCF(_coarseTrialSpace.getSize(), _fineTrialSpace.getSize());
    auto &matrixMgr = DiscreteSpaceMatrixMgr_1D::get();
    LHCF = matrixMgr.getData(_coarseTestSpace, _fineTestSpace);
    RHCF = matrixMgr.getData(_coarseTrialSpace, _fineTrialSpace);
    auto LHCFIndices = matrixMgr.nonZeroInCol(_coarseTestSpace, _fineTestSpace, i);
    auto RHCFIndices = matrixMgr.nonZeroInCol(_coarseTestSpace, _fineTestSpace, j);
    BEM::Complex coarseData{0.0};
    for (const int row : RHCFIndices) {
        for (const int col : LHCFIndices) {
            coarseData += _coarseOperator.indexedDuality(col, row)*LHCF.col(i)[col]*RHCF.col(j)[row];
        }
    }
    return fineData - coarseData;
}

/**
 * @brief: Efficient computation of singular matrix terms for later approximation via MEIM.
 */
template <class OperatorType>
double MultiLevelBEMOperator<OperatorType>::error()
{
    assert(_massMatrix);
    return (_massMatrix->operatorNorm())/_fineOperator.getMatrix().operatorNorm();
}

#endif
