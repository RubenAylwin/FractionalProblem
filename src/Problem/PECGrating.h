#ifndef PEC_GRATING
#define PEC_GRATING

#include <ParametricProblem.h>
#include <Curve.h>
#include <MyTypes.h>
#include <memory>
#include <vector>
#include <random>
#include <string>

class DiscreteSpaceOnCurve_1D;
class Operator;
class PeriodicCurve;

/////////////////////////////
// Classes for PEC problem //
/////////////////////////////

//TODO: Change operators and spaces to spaces over meshes instead of partitions of [0,1]

/**
 * @brief: PEC grating problem
 */
class PECGrating : public QPProblem
{
public:
    PECGrating(PeriodicCurve* curve, double angle, double wavenumber, unsigned numElements, unsigned greenTerms);
    ~PECGrating(void);
    void buildDiscrete(void) override;
    virtual double efficiency(int mode) const;
    virtual double matrixError(void) {return 0.0;};
    virtual const BEM::Matrix &interpolatedMatrix(void) {return *_matrix;};
    const PeriodicCurve &getCurve(void) {return *_curve;};
    const DiscreteSpaceOnCurve_1D &getSpace(void) override {return *_space;};
protected:
    std::unique_ptr<PeriodicCurve> _curve = nullptr;
    std::unique_ptr<DiscreteSpaceOnCurve_1D> _space = nullptr;
    std::unique_ptr<Operator> _V;
    unsigned _numElements = 0U;
};

/**
 * @brief: PEC grating factory. Randomized over trigonometric geometries.
 */
class UQProblemFactory : public RandomProblemFactory
{
public:
    UQProblemFactory(double period, double angle, double wavenumber, double pertSize, int dimDecay, unsigned numElements, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients, unsigned greenTerms = 50U);
    ~UQProblemFactory();
protected:
    std::unique_ptr<Problem> generateNewProblem(unsigned counter) const override;
    std::unique_ptr<Problem> generateNewTestingProblem(unsigned counter) const override;
    double _period;
    double _angle;
    double _wavenumber;
    double _pertSize;
    int _dimDecay;
    unsigned _numElements;
    std::vector<double> _sineCoefficients;
    std::vector<double> _cosineCoefficients;
    unsigned _greenTerms = 50;
    std::unique_ptr<PeriodicCurve> _curve = nullptr;
    mutable std::uniform_real_distribution<double> _unif;
    mutable std::mt19937 _rng;
};

#endif
