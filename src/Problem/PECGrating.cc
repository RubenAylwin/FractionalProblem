#include <PECGrating.h>
#include <RegularP0.h>
#include <WeaklySingular.h>
#include <Curve.h>
#include <GreenQP.h>
#include <cmath>
#include <Utilities.h>
#include <cassert>
#include <iostream>
#include <Identity.h>
#include <EmpiricalInterpolation.h>
#include <DiscreteSpaceMatrixMgr.h>
#include <MultiLevelIntegralOperator.h>
#include <Msg.h>
#include <Halton.h>
#include <Sobol.h>
useMessages("PEC_PROB");

/**
 * @brief: Constructor for holographic uq problem. Specified by problem parameters and vectors describing the geometry.
 */
PECGrating::PECGrating(PeriodicCurve* curve, double angle, double wavenumber, unsigned numElements, unsigned greenTerms) :
    QPProblem(curve->getPeriod(), angle, wavenumber),
    _curve{std::unique_ptr<PeriodicCurve>(curve)},
    _numElements{numElements}
{
    _green->setWindowTerms(greenTerms);
    _space = std::unique_ptr<DiscreteSpaceOnCurve_1D>(new RegularP0_1D(*_curve, _numElements));
    _V.reset(new WeaklySingular(*_space, *_space, *_green));
}

/**
 * @brief: Explicit destructor.
 */
PECGrating::~PECGrating(void)
{
}

/**
 * @brief: Build the matrix and rhs for the discrete problem.
 */
void PECGrating::buildDiscrete(void)
{
    msg(0) << "start::PECGrating::buildDiscrete" << endMsg;
    _V->assembleMassMatrix();
    _matrix.reset(new BEM::Matrix(_V->getMatrix()));

    auto planeWave = BEM::generatePlaneWave(_angle, _wavenumber, *_curve);
    BEM::plotFunction("pw", *planeWave);
    _rhs.reset(new BEM::ColVector(_space->testAgainstBasis(*planeWave)));
    msg(1) << "MatrixDim: " << _matrix->rows() << "x" << _matrix->cols() << endMsg;
    msg(1) << "rhsDim: " << _rhs->size() << endMsg;
    msg(0) << "end::PECGrating::buildDiscrete" << endMsg;
}

/**
 * @brief: Compute the diffraction efficiency for a given mode.
 */
double PECGrating::efficiency(int mode) const
{
    double Bn = _wavenumber*std::sin(_angle) + mode*2.0*M_PI/(_period);
    // If the mode is evanescent, just early return 0.
    if (std::abs(_wavenumber) < std::abs(Bn)) {
        return 0.0;
    }
    
    double Gn = std::sqrt(std::pow(_wavenumber, 2) - std::pow(Bn, 2));

    // Associated plane wave
    auto effPW = BEM::generatePlaneWave(BEM::I*Bn, BEM::I*Gn, *_curve);
    auto vec =  _space->testAgainstBasis(*effPW);

    // Final computation
    BEM::Complex Un = _solution->dot(vec)/(2*Gn);
    return std::pow(std::abs(Un), 2)*Gn/(_period*_period*_wavenumber*std::cos(_angle));
}

/**
 * @brief: Construction for a factory that generates UQ problems with specified random parameters.
 */
UQProblemFactory::UQProblemFactory(double period, double angle, double wavenumber, double pertSize, int dimDecay, unsigned numElements, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients, unsigned greenTerms) :
    RandomProblemFactory(),
    _period{period},
    _angle{angle},
    _wavenumber{wavenumber},
    _pertSize{pertSize},
    _dimDecay{dimDecay},
    _numElements{numElements},
    _sineCoefficients{sineCoefficients},
    _cosineCoefficients{cosineCoefficients},
    _greenTerms{greenTerms},
    _curve{nullptr},
    _unif(-1., 1.),
    _rng((std::random_device())())
{
    _curve.reset(new TrigonometricCurve(period, 0.0, std::vector<double>(_sineCoefficients), std::vector<double>(_cosineCoefficients)));
}

/**
 * @brief: Explicit destructor.
 */
UQProblemFactory::~UQProblemFactory(void)
{
}

//Change to 1 to use halton numbers on the coming function.
#define USE_HALTON 0
/**
 * @brief: Generate a problem.
 * @note: Counter given because implementation has changed between purely random and Halton points.
 */
std::unique_ptr<Problem> UQProblemFactory::generateNewProblem([[maybe_unused]] unsigned counter) const
{
    std::vector<double> randomNumber(_sineCoefficients.size() + _cosineCoefficients.size(), 0);
    # if USE_HALTON
    auto haltonResult = halton(counter, randomNumber.size());
    # endif
    for (size_t i = 0; i < randomNumber.size(); ++i) {
        # if USE_HALTON
        // This part takes points from a halton series.
        randomNumber[i] = haltonResult[i];
        # else
        // This part creates purely random points from a uniform distribution
        randomNumber[i] = _unif(_rng);
        # endif
    }
    # if USE_HALTON
    delete[] haltonResult;
    # endif
    std::vector<double> sineCoefficients(_sineCoefficients.size(), 0);
    std::vector<double> cosineCoefficients(_cosineCoefficients.size(), 0);

    size_t i = 0;
    for (; i < sineCoefficients.size(); ++i) {
        sineCoefficients[i] = _sineCoefficients[i] + (randomNumber[i])*_pertSize/(std::pow(i+1., _dimDecay));
    }
    for (; i < sineCoefficients.size() + cosineCoefficients.size() ; ++i) {
        size_t j = i - sineCoefficients.size();
        cosineCoefficients[j] = _cosineCoefficients[j] + (randomNumber[i])*_pertSize/std::pow(j+1., _dimDecay);
    }
    
    auto ptr = new PECGrating(new TrigonometricCurve(_period, 0., std::move(sineCoefficients), std::move(cosineCoefficients)), _angle, _wavenumber, _numElements, _greenTerms);
    return std::unique_ptr<Problem>(ptr);
}



//Change to 1 to use sobol numbers on the coming function.
#define USE_SOBOL 0
/**
 * @brief: Generate a problem for testing.
 * @desc: Essentialy the same as before. Was made to be able to generate points from a different origin (either random again or Sobol instead of Halton).
 */
std::unique_ptr<Problem> UQProblemFactory::generateNewTestingProblem([[maybe_unused]] unsigned counter) const
{
    std::vector<double> sobolNumber(_sineCoefficients.size() + _cosineCoefficients.size(), 0);
    # if USE_SOBOL
    double *sobolResult = new double[sobolNumber.size()];
    auto newCounter = static_cast<long long int>(counter);
    i8_sobol(sobolNumber.size(), &newCounter, sobolResult);
    # endif
    for (size_t i = 0; i < sobolNumber.size(); ++i) {
        # if USE_SOBOL
        sobolNumber[i] = 1. - 2.*sobolResult[i];
        # else
        sobolNumber[i] = _unif(_rng);
        # endif
    }
    #if USE_SOBOL
    delete[] sobolResult;
    # endif
    std::vector<double> sineCoefficients(_sineCoefficients.size(), 0);
    std::vector<double> cosineCoefficients(_cosineCoefficients.size(), 0);

    size_t i = 0;
    for (; i < sineCoefficients.size(); ++i) {
        sineCoefficients[i] = _sineCoefficients[i] + (sobolNumber[i])*_pertSize/(std::pow(i+1., _dimDecay));
    }
    for (; i < sineCoefficients.size() + cosineCoefficients.size() ; ++i) {
        int j = i - sineCoefficients.size();
        cosineCoefficients[j] = _cosineCoefficients[j] + (sobolNumber[i])*_pertSize/std::pow(j+1., _dimDecay);
    }
    
    auto ptr = new PECGrating(new TrigonometricCurve(_period, 0., std::move(sineCoefficients), std::move(cosineCoefficients)), _angle, _wavenumber, _numElements, _greenTerms);
    return std::unique_ptr<Problem>(ptr);
}
