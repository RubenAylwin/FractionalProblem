#include <L2Operator.h>
#include <DiscreteSpace.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>

L2::L2(ExplicitScalarFunction_1D &function, const DiscreteSpaceOnCurve_1D &trialSpace, DiscreteSpaceOnCurve_1D &testSpace) :
    DifferentialOperator(trialSpace, testSpace),
    _kernel{function}

{
    _integrator1D.reset(new GaussLegendre_1D(10));
}

void L2::assembleMassMatrix(void)
{
    if (_massMatrix) {
        return;
    }
    _massMatrix.reset(new BEM::Matrix(BEM::Matrix::Zero(Operator::_testSize, Operator::_trialSize)));
    auto &testElements = _testSpace.getElements();
    auto &trialElements = _testSpace.getElements();
    assert(testElements.size()==trialElements.size());
    auto &partition = _testSpace.getPartition();
    for (size_t element = 0; element < testElements.size(); ++element) {
        auto &testElement = testElements[element];
        auto &trialElement = trialElements[element];
        assert(trialElement == testElement);
        for (unsigned i_test = 0; i_test < testElement.getDof(); ++i_test) {
            for (unsigned i_trial = 0; i_trial < trialElement.getDof(); ++i_trial) {
                auto testIndex = _testSpace.globalIndex(element, i_test);
                auto trialIndex = _trialSpace.globalIndex(element, i_trial);
                auto &test = _testSpace.basisFunction(testIndex);
                auto &trial = _trialSpace.basisFunction(trialIndex);
                auto product = test*trial;
                (*_massMatrix)(testIndex, trialIndex) += _integrator1D->integrate(partition[testElement.getA()], partition[testElement.getB()], (*product*_kernel));
            }
        }
    }    
}

BEM::Complex L2::indexedDuality(const unsigned i, const unsigned j)
{
    auto &test = _testSpace.basisFunction(i);
    auto &trial = _trialSpace.basisFunction(j);
    auto product = test*trial;
    return  _integrator1D->integrate(test.support().first, test.support().second, (*product*_kernel));
}
