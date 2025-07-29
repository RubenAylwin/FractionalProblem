#include <PECGratingReduced.h>
#include <EmpiricalInterpolation.h>
#include <Utilities.h>
#include <iostream>
#include <progressbar.hpp>
#include <DiscreteSpaceMatrixMgr.h>
#include <DiscreteSpace.h>
#include <GreenQP.h>
#include <IntegralOperator.h>
#include <Msg.h>

useMessages("UQ_RB_PROB");


/**
 * @brief: Problem constructor. Requires a reference to a trained empirical interpolator.
 */
PECGratingReduced::PECGratingReduced(PeriodicCurve* curve, double angle, double wavenumber, unsigned numElements, const BEM::EmpiricalInterpolationPtr &empiricalGreen, const BEM::EmpiricalInterpolationPtr &empiricalRhs,  const std::shared_ptr<std::vector<BEM::Matrix>> &greenMatrices,  const std::shared_ptr<std::vector<BEM::ColVector>> &rhsVectors , const std::shared_ptr<BEM::Matrix> &bemPodBasis, unsigned greenTerms) :
    PECGrating(curve, angle, wavenumber, numElements, greenTerms),
    _empiricalGreen{empiricalGreen},
    _empiricalRhs{empiricalRhs},
    _greenMatrices{greenMatrices},
    _rhsVectors{rhsVectors},
    _bemPodBasis{bemPodBasis}
{
}

/**
 * @brief: Explicit destructor
 */
PECGratingReduced::~PECGratingReduced(void)
{
}

/**
 * @brief: Build a reduced problem.
 */
void PECGratingReduced::buildDiscrete(void)
{
    msg(5) << "(start)PECGratingReduced::buildDiscrete"<< endMsg;
    {
        auto indices = _empiricalGreen->getInterpolationIndices();
        std::vector<BEM::Complex> data{};
        for (size_t i = 0; i < indices.size(); ++i) {
            const auto &index = indices[i];
            int col = (int)(index/_numElements);
            int row = index - col*_numElements;
            data.push_back(_V->indexedDuality(row, col));
        }

        auto solution = _empiricalGreen->getSolution(data);
        auto &greenMatrices = *_greenMatrices;
        BEM::Matrix matrix = greenMatrices[0]*solution[0];
        for (size_t i = 1; i < indices.size(); ++i) {
            matrix += greenMatrices[i]*solution[i];
        }
        _matrix.reset(new BEM::Matrix(std::move(matrix)));
    }

    {
        auto planeWave = BEM::generatePlaneWave(_angle, _wavenumber, *_curve);
        auto indices = _empiricalRhs->getInterpolationIndices();
        std::vector<BEM::Complex> data{};
        for (size_t i = 0; i < indices.size(); ++i) {
            const auto &index = indices[i];
            data.push_back(_space->testAgainstBaseElement(*planeWave, index));
        }

        auto solution = _empiricalRhs->getSolution(data);
        auto &rhsVectors = *_rhsVectors;
        BEM::ColVector rhs = rhsVectors[0]*solution[0];
        for (size_t i = 1; i < indices.size(); ++i) {
            rhs += rhsVectors[i]*solution[i];
        }
        _rhs.reset(new BEM::ColVector(std::move(rhs)));
    }
    msg(6) << "MatrixDim: " << _matrix->rows() << "x" << _matrix->cols() << endMsg;
    msg(6) << "RhsDim: " << _rhs->size() << endMsg;
    msg(5) << "(end)PECGratingReduced::buildDiscrete"<< endMsg;
}


int PECGratingReduced::solve(void)
{
    msg(5) << "(start)PECGratingReduced::solve" << endMsg;
    auto first = BEM::now();
    BEM::ColVector podSolution = _matrix->colPivHouseholderQr().solve(*_rhs);
    auto second = BEM::now();
    msg(6) << "Solved linear system" << endMsg;

    BEM::ColVector solution = podSolution[0]*_bemPodBasis->col(0);
    msg(6) << "Initialized solution" << endMsg;
    for (int i = 1; i < podSolution.size(); ++i) {
        msg(6) << "i "  << podSolution.size() << endMsg;
        solution += podSolution[i]*_bemPodBasis->col(i);
    }
    msg(6) << "Finished solution" << endMsg;
    _solution.reset(new BEM::ColVector(solution));
    msg(5) << "(end)PECGratingReduced::solve" << endMsg;
    return BEM::timeDifference(first, second);
}

/**
 * @brief: Check the error between the real matrix and its interpolation.
 */
double PECGratingReduced::matrixError(void)
{
    _V->assembleMassMatrix();
    
    auto indices = _empiricalGreen->getInterpolationIndices();
    std::vector<BEM::Complex> data{};
    for (size_t i = 0; i < indices.size(); ++i) {
        const auto &index = indices[i];
        int col = (int)(index/_numElements);
        int row = index - col*_numElements;
        data.push_back(_V->indexedDuality(row, col));
    }
    auto interpolation = (BEM::compress(_empiricalGreen->interpolate(data)));
    
    return (_V->getMatrix() - interpolation).operatorNorm();
}

/**
 * @brief: Get the interpolated matrix.
 */
const BEM::Matrix &PECGratingReduced::interpolatedMatrix(void)
{
    auto indices = _empiricalGreen->getInterpolationIndices();
    std::vector<BEM::Complex> data{};
    for (size_t i = 0; i < indices.size(); ++i) {
        const auto &index = indices[i];
        int col = (int)(index/_numElements);
        int row = index - col*_numElements;
        data.push_back(_V->indexedDuality(row, col));
    }
    _matrix.reset(new BEM::Matrix(BEM::compress(_empiricalGreen->interpolate(data))));
    return *_matrix;
}

/**
 * @brief: Get the RHS.
 */
const BEM::ColVector &PECGratingReduced::getRHS(void)
{
    auto indices = _empiricalRhs->getInterpolationIndices();
    auto planeWave = BEM::generatePlaneWave(_angle, _wavenumber, *_curve);
    _rhs.reset(new BEM::ColVector(_space->testAgainstBasis(*planeWave)));
    return *_rhs;
    std::vector<BEM::Complex> data{};
    for (size_t i = 0; i < indices.size(); ++i) {
        const auto &index = indices[i];
        data.push_back(_space->testAgainstBaseElement(*planeWave, index));
    }
    _rhs.reset(new BEM::ColVector(_empiricalRhs->interpolate(data)));
    return *_rhs;
}

/**
 * @brief: Get error between interpolated and real RHS.
 */
double PECGratingReduced::rhsError(void)
{

    auto indices = _empiricalRhs->getInterpolationIndices();
    auto planeWave = BEM::generatePlaneWave(_angle, _wavenumber, *_curve);
        
    std::vector<BEM::Complex> data{};
    for (size_t i = 0; i < indices.size(); ++i) {
        const auto &index = indices[i];
        data.push_back(_space->testAgainstBaseElement(*planeWave, index));
    }
    auto interpolation = (_empiricalRhs->interpolate(data));
    auto realVector = _space->testAgainstBasis(*planeWave);
    return (interpolation - realVector).norm();
}

/**
 * @brief: Get error between real and RB solution.
 */
double PECGratingReduced::solutionError(void)
{
    msg(5) << "(start)PECGratingReduced::solutionError" << endMsg;
    buildDiscrete();
    solve();
    auto& basis = *_bemPodBasis;
    auto& rbSolution = *_solution;
    msg(6) << "Bases: " << basis.cols() << endMsg;
    msg(6) << "SolLen: " << rbSolution.size() << endMsg;
    BEM::ColVector expandedSolution = basis.col(0)*rbSolution[0];
    for (int i = 1; i < rbSolution.size(); ++i) {
        expandedSolution += basis.col(i)*rbSolution[i];
    }

    _V->assembleMassMatrix();
    auto planeWave = BEM::generatePlaneWave(_angle, _wavenumber, *_curve);
    auto rhs = BEM::ColVector(_space->testAgainstBasis(*planeWave));
    BEM::ColVector solution = _V->getMatrix().colPivHouseholderQr().solve(rhs);

    msg(5) << "(end)PECGratingReduced::solutionError" << endMsg;
    return (solution - expandedSolution).norm()/solution.norm();
}

/**
 * @brief: Constructor.
 */
UQProblemFactoryReduced::UQProblemFactoryReduced(double period, double angle, double wavenumber, double pertSize, int dimDecay, unsigned numElements, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients, unsigned greenTerms) :
    UQProblemFactory(period, angle, wavenumber, pertSize, dimDecay, numElements, sineCoefficients, cosineCoefficients, greenTerms)
{
}

/**
 * @brief: Explicit destructor.
 */
UQProblemFactoryReduced::~UQProblemFactoryReduced(void)
{
}

/**
 * @brief: Initialize necessary structures for SVD.
 */
void UQProblemFactoryReduced::initializeSVD(void)
{
    msg(5) << "(start)UQProblemFactoryReduced::initializeSVD" << endMsg;
    _trained = false;
    _greenMatrices.reset(new std::vector<BEM::Matrix>);
    _rhsVectors.reset(new std::vector<BEM::ColVector>);
    _empiricalGreen.reset(new EmpiricalInterpolation("Green", _numElements*_numElements));
    _empiricalRhs.reset(new EmpiricalInterpolation("RHS", _numElements));
    _empiricalSolution.reset(new EmpiricalInterpolation("Solution", _numElements));
    msg(5) << "(end)UQProblemFactoryReduced::initializeSVD" << endMsg;
    
}

/**
 * @brief: Sample the space with MC.
 * @desc: Now, depending on some macro variables, generateNewProblem may use Halton points.
 */
void UQProblemFactoryReduced::sampleMC(unsigned samples)
{
    msg(5) << "(start)UQProblemFactoryReduced::sampleMC" << endMsg;
    initializeSVD();
    progressbar bar(samples);
    for (unsigned sample = 0; sample < samples; ++sample) {
        bar.update();
        //Generate and solve a problem.
        auto prob = RandomProblemFactory::generateNewProblem();
        prob->buildDiscrete();
        prob->solve();
        //Save the matrix, rhs and solution to empirical interpolators.
        _empiricalGreen->loadTrainingData(BEM::stretch(prob->getMatrix()));
        _empiricalRhs->loadTrainingData(BEM::ColVector(prob->getRHS()));
        _empiricalSolution->loadTrainingData(BEM::ColVector(prob->getSolutionVec()));
    }
    std::clog << std::endl;
    msg(5) << "(end)UQProblemFactoryReduced::sampleMC" << endMsg;
}

/**
 * @brief: Take the data from the interpolators and save it to the factory to have the necessary bases.
 */
void UQProblemFactoryReduced::saveEIData(void)
{
    msg(5) << "(start)UQProblemFactoryReduced::saveEIData" << endMsg;
    _bemPodBasis.reset(new BEM::Matrix(_empiricalSolution->podBasis()));
    
    const auto &greenPoD = _empiricalGreen->podBasis();
    for (int col = 0; col < greenPoD.cols(); ++col) {
        _greenMatrices->push_back(BEM::product(BEM::product(_bemPodBasis->adjoint(), BEM::compress(greenPoD.col(col))), *_bemPodBasis));
    }

    const auto &rhsPod = _empiricalRhs->podBasis();
    for (int col = 0; col < rhsPod.cols(); ++col) {
        _rhsVectors->push_back(BEM::product(_bemPodBasis->adjoint(), rhsPod.col(col)));
    }
    _trained = true;
    msg(5) << "(end)UQProblemFactoryReduced::saveEIData" << endMsg;
}

/**
 * @brief: Train the models to given tolerance.
 */
void UQProblemFactoryReduced::trainEmpirical(unsigned samples, double tolerance)
{
    sampleMC(samples);
    //Matrix and RHS get lower tolerance so that their error does not overly influence the solution error.
    _empiricalGreen->buildInterpolator(tolerance/10./*Tolerance*/);
    msg(2) << "Built PoD-Green interpolator" << endMsg;

    _empiricalRhs->buildInterpolator(tolerance/10./*Tolerance*/);
    msg(2) << "Built PoD-Rhs interpolator" << endMsg;

    _empiricalSolution->buildInterpolator(tolerance/*Tolerance*/);
    msg(2) << "Built PoD-Solution interpolator" << endMsg;

    saveEIData();
}

/**
 * @brief: Train until number of bases.
 */
void UQProblemFactoryReduced::trainEmpirical(unsigned samples, unsigned bases)
{
    sampleMC(samples);
    if (bases > samples) {
        bases = samples;
    }
    if (bases > _numElements*_numElements) {
        bases = _numElements*_numElements;
    }
    _empiricalGreen->buildInterpolator(bases);
    msg(2) << "Built PoD-Green interpolator" << endMsg;

    // The green function usually requires more bases than the RHS and solution.
    // Sometimes more than the whole discrete dimension N.
    if (bases > _numElements) {
        bases = _numElements;
    }

    _empiricalRhs->buildInterpolator(bases);
    msg(2) << "Built PoD-Rhs interpolator" << endMsg;

    _empiricalSolution->buildInterpolator(bases);
    msg(2) << "Built PoD-Solution interpolator" << endMsg;
    
    saveEIData();
}

/**
 * @brief: Generate a new problem.
 */
std::unique_ptr<Problem> UQProblemFactoryReduced::generateNewProblem(unsigned counter) const
{
    // If not trained, get a real problem from the base factory.
    if (not _trained) {
        return UQProblemFactory::generateNewProblem(counter);
    }

    // If trained, get a reduced problem.
    std::vector<double> mcData(_sineCoefficients.size() + _cosineCoefficients.size(), 0);
    for (size_t i = 0; i < mcData.size(); ++i) {
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
    ptr = new PECGratingReduced(curve, _angle, _wavenumber, _numElements, _empiricalGreen, _empiricalRhs, _greenMatrices, _rhsVectors, _bemPodBasis,  _greenTerms);

    return std::unique_ptr<Problem>(ptr);
}

/**
 * @brief: Generate a (reduced) problem from a given curve (not the parameters).
 */
std::unique_ptr<PECGrating> UQProblemFactoryReduced::generateProblemFromCurve(PeriodicCurve *curve)
{
    assert(_trained);
    PECGrating *ptr = new PECGratingReduced(curve, _angle, _wavenumber, _numElements, _empiricalGreen, _empiricalRhs, _greenMatrices, _rhsVectors, _bemPodBasis,  _greenTerms);
    return std::unique_ptr<PECGrating>(ptr);
}
