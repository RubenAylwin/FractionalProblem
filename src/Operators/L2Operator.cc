#include <L2Operator.h>
#include <Utilities.h>
#include <DiscreteSpace.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <ScalarValuedFunction.h>
#include <Mesh.h>

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


L2Mesh::L2Mesh(const DiscreteSpaceMesh &trial, const DiscreteSpaceMesh &test):
    DifferentialOperatorMesh(trial, test)
{
    _integrator1D.reset(new GaussLegendre_1D(5));//Hardcoded. Should be enough (maybe even 2 in enough).
}

BEM::Complex L2Mesh::indexedDuality(const unsigned i, const unsigned j)
{
    auto &test = _testSpace.basisFunction(i);
    auto &trial = _trialSpace.basisFunction(j);
    auto &mesh = _testSpace.getMesh();
    BEM::Complex result = 0.0;
    unsigned l = 0;    
    auto integrationElementSize = mesh.getElement(0).getSize();
    ExplicitScalarFunction_1D integrand([&](double t) {
        auto point = mesh.getElement(l)(t);
        return test.evaluate(l, t)*trial.evaluate(l, t)*integrationElementSize*(*_dFun)(point.getX());
    }
        );
    // Right now we integrate over all elements, but we know that the integration domain is less that this.
    for (unsigned k : trial.support()) {
        l=k;
        integrationElementSize = mesh.getElement(l).getSize();
        result += _integrator1D->integrate(0, 1, integrand);
    }
    return result;
}
    
