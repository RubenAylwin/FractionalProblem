#include <PECGratingReducedML.h>
#include <RegularP0.h>
#include <EmpiricalInterpolation.h>
#include <Utilities.h>
#include <iostream>
#include <progressbar.hpp>
#include <DiscreteSpaceMatrixMgr.h>
#include <DiscreteSpace.h>
#include <GreenQP.h>
#include <IntegralOperator.h>
#include <Msg.h>
#include <WeaklySingular.h>
#include <MultiLevelIntegralOperator.h>

useMessages("UQ_RB_ML_PROB");

/////////////////////////////////////////////////////////////////////
// Under development.                                              //
// Idea is to have faster training by training at different levels //
/////////////////////////////////////////////////////////////////////

PECGratingML::PECGratingML(PeriodicCurve* curve, double angle, double wavenumber, unsigned numElementsHigh, unsigned numElementsLow, unsigned greenTerms) :
    PECGrating(curve, angle, wavenumber, numElementsHigh, greenTerms),
    _numElementsLow(numElementsLow)
{
    if (_numElementsLow > 0) {
        _coarseSpace.reset(new RegularP0_1D(*_curve, _numElementsLow));
        _V.reset(new MultiLevelBEMOperator<WeaklySingular>(*_space, *_space, *_coarseSpace, *_coarseSpace, *_green));
    }
}

PECGratingML::~PECGratingML(void)
{
}

void PECGratingML::buildDiscrete(void)
{
    PECGrating::buildDiscrete();
    if (_numElementsLow > 0) {
        auto& matrixMgr = DiscreteSpaceMatrixMgr_1D::get();
        auto LHFC = matrixMgr.getData(*_coarseSpace, *_space).transpose();
        auto planeWave = BEM::generatePlaneWave(_angle, _wavenumber, *_curve);
        _coarseRhs.reset(new BEM::ColVector(_coarseSpace->testAgainstBasis(*planeWave)));
        _fineRhs.reset(new BEM::ColVector(*_rhs));
        _rhs.reset(new BEM::ColVector(*_fineRhs - LHFC*(*_coarseRhs)));
    }
}

int PECGratingML::solve(void)
{
    if (_numElementsLow == 0) {
        return PECGrating::solve();
    }
    auto first = BEM::now();
    auto& matrixMgr = DiscreteSpaceMatrixMgr_1D::get();
    auto LHFC = matrixMgr.getData(*_coarseSpace, *_space).transpose();
    auto V = static_cast<MultiLevelBEMOperator<WeaklySingular>*>(_V.get());
    _solution.reset(new BEM::ColVector(V->getFineMatrix().colPivHouseholderQr().solve(*_fineRhs) - LHFC*(V->getCoarseMatrix().colPivHouseholderQr().solve(*_coarseRhs))));
    auto second = BEM::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(second - first).count();
}


PECGratingReducedML::PECGratingReducedML(PeriodicCurve* curve, double angle, double wavenumber, std::vector<unsigned> numElementsVector, std::shared_ptr<std::vector<EmpiricalInterpolation>> empiricalGreen, std::shared_ptr<std::vector<EmpiricalInterpolation>> empiricalRhs, std::shared_ptr<BEM::Matrix> bemPodMatrix, std::shared_ptr<std::vector<BEM::Matrix>> bemMatrices, std::shared_ptr<std::vector<BEM::ColVector>> bemRhs, unsigned greenTerms) :
    PECGrating(curve, angle, wavenumber, numElementsVector.back(), greenTerms),
    _numElementsVector(numElementsVector),
    _empiricalGreen{empiricalGreen},
    _empiricalRhs{empiricalRhs},
    _bemPodMatrix{bemPodMatrix},
    _bemMatrices{bemMatrices},
    _bemRhs{bemRhs},
    _spaces{},
    _operators{}
{
    assert(_numElementsVector[0] == 0);
    for (unsigned i = 1; i < _numElementsVector.size(); ++i) {
        _spaces.emplace_back(new RegularP0_1D(*_curve, _numElementsVector[i]));
    }
    _operators.emplace_back(new WeaklySingular(*(_spaces[0]), *(_spaces[0]), *_green));
    for (unsigned i = 1; i < _spaces.size(); ++i) {
        _operators.emplace_back(new MultiLevelBEMOperator<WeaklySingular>(*(_spaces[i]), *(_spaces[i]), *(_spaces[i-1]), *(_spaces[i-1]), *_green));
    }
}

PECGratingReducedML::~PECGratingReducedML(void)
{
}

void PECGratingReducedML::interpolateGreen(void)
{
    auto &matrixMgr = DiscreteSpaceMatrixMgr_1D::get();
    auto &fine = *(_spaces.back());
    for (unsigned level = 0; level < _operators.size(); ++level) {
        auto LHFC = matrixMgr.getData(*_spaces[level], fine).transpose();
        auto RHFC = matrixMgr.getData(*_spaces[level], fine);
        auto &empiricalGreen = _empiricalGreen->at(level);
        auto &op = *(_operators[level]);
        auto indices = empiricalGreen.getInterpolationIndices();
        std::vector<BEM::Complex> data;
        for (unsigned index = 0; index < indices.size(); ++index) {
            int col = (int)(indices[index]/_numElementsVector[level+1]);
            int row = indices[index] - col*_numElementsVector[level+1];
            data.push_back(op.indexedDuality(row, col));
        }
        if (not _matrix) {
            //Coarse level
            _matrix.reset(new BEM::Matrix(BEM::product(BEM::product(LHFC, BEM::compress(empiricalGreen.interpolate(data))), RHFC)));
            continue;
        }
        *_matrix += BEM::product(BEM::product(LHFC, BEM::compress(empiricalGreen.interpolate(data))), RHFC);
    }
}

double PECGratingReducedML::matrixError(void)
{
    interpolateGreen();
    auto &fullMat = static_cast<MultiLevelBEMOperator<WeaklySingular>*>(_operators.back().get())->getFineMatrix();
    return (*_matrix - fullMat).operatorNorm();
}

const BEM::Matrix &PECGratingReducedML::interpolatedMatrix(void)
{
    interpolateGreen();
    return *_matrix;
}

UQProblemFactoryReducedML::UQProblemFactoryReducedML(double period, double angle, double wavenumber, double pertSize, int dimDecay, std::vector<unsigned> numElementsVector, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients, unsigned greenTerms):
    UQProblemFactoryReduced(period, angle, wavenumber, pertSize, dimDecay, numElementsVector.back(), sineCoefficients, cosineCoefficients, greenTerms),
    _numElementsVector{numElementsVector},
    _empiricalGreen{std::make_shared<std::vector<EmpiricalInterpolation>>()},
    _empiricalRhs{std::make_shared<std::vector<EmpiricalInterpolation>>()},
    _empiricalSolution{std::make_shared<std::vector<EmpiricalInterpolation>>()},
    _spaces()
{
    msg(5) << "(start)UQProblemFactoryReducedML::UQProblemFactoryReducedML" << endMsg;
    if (_numElementsVector[0] != 0) {
        _numElementsVector.insert(_numElementsVector.begin(), 0U);
    }
    assert(_numElementsVector[0] == 0);
    for (unsigned i = 1; i < _numElementsVector.size(); ++i) {
        _spaces.emplace_back(new RegularP0_1D(*_curve, _numElementsVector[i]));
    }
    msg(5) << "(end)UQProblemFactoryReducedML::UQProblemFactoryReducedML" << endMsg;
}


UQProblemFactoryReducedML::~UQProblemFactoryReducedML(void)
{
}

void UQProblemFactoryReducedML::initializeSVD(void)
{
    msg(5) << "(start)UQProblemFactoryReducedML::initializeSVD"<< endMsg;
    _trained = false;
    for (unsigned i = 1; i < _numElementsVector.size(); ++i) {
        _empiricalGreen->emplace_back((_numElementsVector[i]*_numElementsVector[i]));
        _empiricalRhs->emplace_back((_numElementsVector[i]));
        _empiricalSolution->emplace_back((_numElementsVector[i]));
    }
    msg(5) << "(end)UQProblemFactoryReducedML::initializeSVD" << endMsg;
}

PECGratingML UQProblemFactoryReducedML::generateNewTrainingProblem(unsigned level, [[maybe_unused]] unsigned counter)
{
    std::vector<double> mcData(_sineCoefficients.size() + _cosineCoefficients.size(), 0);
    for (unsigned i = 0; i < mcData.size(); ++i) {
        mcData[i] = _unif(_rng);
    }

    std::vector<double> sineCoefficients(_sineCoefficients.size(), 0);
    std::vector<double> cosineCoefficients(_cosineCoefficients.size(), 0);

    size_t i = 0;
    for (; i < sineCoefficients.size(); ++i) {
        sineCoefficients[i] = _sineCoefficients[i] + (1.0 - 2.*mcData[i])*_pertSize/(std::pow(i+1., _dimDecay));
    }
    for (; i < sineCoefficients.size() + cosineCoefficients.size() ; ++i) {
        int j = i - sineCoefficients.size();
        cosineCoefficients[j] = _cosineCoefficients[j] + (1.0 - 2.*mcData[i])*_pertSize/std::pow(j+1., _dimDecay);
    }

    auto curve = new TrigonometricCurve(_period, 0., std::move(sineCoefficients), std::move(cosineCoefficients));
    return PECGratingML(curve, _angle, _wavenumber, _numElementsVector[level+1], _numElementsVector[level], _greenTerms);
}

void UQProblemFactoryReducedML::sampleMC(const std::vector<unsigned> &samples)
{
    msg(5) << "(start)UQProblemFactoryReducedML::sampleMC"<< endMsg;
    for (unsigned level = 0; level < _numElementsVector.size() - 1 ; ++level) {
        for (unsigned sample = 0; sample < samples[level]; ++sample) {
            auto problem = generateNewTrainingProblem(level, sample);
            problem.buildDiscrete();
            problem.solve();
            _empiricalGreen->at(level).loadTrainingData(BEM::stretch(problem.getMatrix()));
            _empiricalRhs->at(level).loadTrainingData(BEM::ColVector(problem.getRHS()));
            _empiricalSolution->at(level).loadTrainingData(BEM::ColVector(problem.getSolutionVec()));
        }
        _empiricalGreen->at(level).computeSVD();
        _empiricalRhs->at(level).computeSVD();
        _empiricalSolution->at(level).computeSVD();
    }
    msg(5) << "(end)UQProblemFactoryReducedML::sampleMC"<< endMsg;
}

void UQProblemFactoryReducedML::trainEmpirical(const std::vector<unsigned> &samples, unsigned bases)
{
    msg(5) << "(start)UQProblemFactoryReducedML::trainEmpirical"<< endMsg;
    initializeSVD();
    assert(samples.size() == _numElementsVector.size() - 1);
    sampleMC(samples);
    unsigned totalDim = 0;
    for (const auto val : _numElementsVector) {
        totalDim+=val;
    }
        
    BEM::getBasesPerLevel(*_empiricalGreen, bases);
    if (bases > totalDim) {
        bases = totalDim;
    }

    BEM::getBasesPerLevel(*_empiricalRhs, bases);
    BEM::getBasesPerLevel(*_empiricalSolution, bases);
    _bemPodMatrix.reset(new BEM::Matrix(_numElementsVector.back(), bases));
    _bemPodMatrix->setZero();

    saveEIData();
    msg(5) << "(end)UQProblemFactoryReducedML::trainEmpirical"<< endMsg;
}

void UQProblemFactoryReducedML::saveEIData(void)
{
    msg(5) << "(start)UQProblemFactoryReducedML::saveEIData" << endMsg;
    auto& matrixMgr = DiscreteSpaceMatrixMgr_1D::get();

    _bemPodBasis.reset(new std::vector<BEM::Matrix>());
    for (unsigned level = 0; level < _numElementsVector.size() - 1; ++level) {
        auto LHFC = matrixMgr.getData(*_spaces[level], *(_spaces.back())).transpose();
        _bemPodBasis->emplace_back(BEM::product(LHFC,_empiricalSolution->at(level).podBasis()));
        msg(6) << "Saved bem matrix pod size " << _bemPodBasis->back().rows() << "x" << _bemPodBasis->back().cols() << endMsg;
    }

    //Fill the PodMatrix
    int col = 0;
    for (unsigned m = 0; m < _bemPodBasis->size(); ++m) {
        for (unsigned j = 0; j < _bemPodBasis->at(m).cols(); ++j) {
            _bemPodMatrix->col(col) = _bemPodBasis->at(m).col(j);
            ++col;
        }
    }
    
    //Fill the green matrices
    _bemMatrices.reset(new std::vector<BEM::Matrix>());
    _bemRedMatrices.reset(new std::vector<BEM::Matrix>());
    for (unsigned level = 0; level < _numElementsVector.size() - 1; ++level) {
        auto LHFC = matrixMgr.getData(*_spaces[level], *(_spaces.back())).transpose();
        auto RHFC = matrixMgr.getData(*_spaces[level], *(_spaces.back()));
        auto &basisLevel = _empiricalGreen->at(level).podBasis();
        for (unsigned col = 0; col < basisLevel.cols(); ++col) {
            _bemMatrices->emplace_back(BEM::product(BEM::product(LHFC, BEM::compress(basisLevel.col(col))), RHFC));
            _bemRedMatrices->emplace_back(BEM::product(BEM::product(_bemPodMatrix->transpose(), _bemMatrices->back()), *_bemPodMatrix));
            msg(6) << "Saved bem matrix size " << _bemMatrices->back().rows() << "x" << _bemMatrices->back().cols() << endMsg;
            msg(6) << "Saved bem reduced matrix size " << _bemRedMatrices->back().rows() << "x" << _bemRedMatrices->back().cols() << endMsg;
        }
    }

    _bemRhs.reset(new std::vector<BEM::ColVector>());
    _bemRedRhs.reset(new std::vector<BEM::ColVector>());
    for (unsigned level = 0; level < _numElementsVector.size() - 1; ++level) {
        auto LHFC = matrixMgr.getData(*_spaces[level], *(_spaces.back())).transpose();
        auto RHFC = matrixMgr.getData(*_spaces[level], *(_spaces.back()));
        auto &basisLevel = _empiricalRhs->at(level).podBasis();
        for (unsigned col = 0; col < basisLevel.cols(); ++col) {
            _bemRhs->emplace_back(LHFC*basisLevel.col(col));
            _bemRedRhs->emplace_back((_bemPodMatrix->transpose()*_bemRhs->back()));
            msg(6) << "Saved bem rhs size " << _bemRhs->back().rows() << "x" << _bemRhs->back().cols() << endMsg;
            msg(6) << "Saved bem reduced rhs size " << _bemRedRhs->back().rows() << "x" << _bemRedRhs->back().cols() << endMsg;
        }
    }
    _trained = true;
    msg(5) << "(end)UQProblemFactoryReducedML::saveEIData" << endMsg;
}

std::unique_ptr<Problem> UQProblemFactoryReducedML::generateNewProblem(unsigned counter) const
{
    if (not _trained) {
        return UQProblemFactory::generateNewProblem(counter);
    }
    
    std::vector<double> mcData(_sineCoefficients.size() + _cosineCoefficients.size(), 0);
    for (unsigned i = 0; i < mcData.size(); ++i) {
        mcData[i] = _unif(_rng);
    }

    std::vector<double> sineCoefficients(_sineCoefficients.size(), 0);
    std::vector<double> cosineCoefficients(_cosineCoefficients.size(), 0);

    size_t i = 0;
    for (; i < sineCoefficients.size(); ++i) {
        sineCoefficients[i] = _sineCoefficients[i] + (1.0 - 2.*mcData[i])*_pertSize/(std::pow(i+1., _dimDecay));
    }
    for (; i < sineCoefficients.size() + cosineCoefficients.size() ; ++i) {
        int j = i - sineCoefficients.size();
        cosineCoefficients[j] = _cosineCoefficients[j] + (1.0 - 2.*mcData[i])*_pertSize/std::pow(j+1., _dimDecay);
    }

    auto curve = new TrigonometricCurve(_period, 0., std::move(sineCoefficients), std::move(cosineCoefficients));
    
    PECGrating *ptr;
    ptr = new PECGratingReducedML(curve, _angle, _wavenumber, _numElementsVector, _empiricalGreen, _empiricalRhs, _bemPodMatrix, _greenMatrices, _rhsVectors,  _greenTerms);

    return std::unique_ptr<Problem>(ptr);

}

std::unique_ptr<PECGrating> UQProblemFactoryReducedML::generateProblemFromCurve(PeriodicCurve *curve)
{
    assert(_trained);
    PECGrating *ptr = new PECGratingReducedML(curve, _angle, _wavenumber, _numElementsVector, _empiricalGreen, _empiricalRhs, _bemPodMatrix, _greenMatrices, _rhsVectors,  _greenTerms);
    return std::unique_ptr<PECGrating>(ptr);
}
