#ifndef PEC_GRATING_REDUCED_ML
#define PEC_GRATING_REDUCED_ML

#include <PECGratingReduced.h>
#include <MyTypes.h>
#include <memory>
#include <vector>
#include <random>

/////////////////////////////////////////////////////////////////////
// Under development.                                              //
// Idea is to have faster training by training at different levels //
/////////////////////////////////////////////////////////////////////

class DiscreteSpaceOnCurve_1D;
class Operator;

class PECGratingML : public PECGrating
{
public:
    PECGratingML(PeriodicCurve* curve, double angle, double wavenumber, unsigned numElementsHigh, unsigned numElementsLow, unsigned greenTerms = 50U);
    ~PECGratingML(void);
    void buildDiscrete(void) override;
    int solve(void) override;
    BEM::ColVector getFineRhs(void) {return *_fineRhs;}
    BEM::ColVector getCoarseRhs(void) {return *_coarseRhs;}
private:
    unsigned _numElementsLow = 0U;
    std::unique_ptr<DiscreteSpaceOnCurve_1D> _coarseSpace = nullptr;
    std::unique_ptr<BEM::ColVector> _fineRhs = nullptr;
    std::unique_ptr<BEM::ColVector> _coarseRhs = nullptr;
};

class PECGratingReducedML : public PECGrating
{
public:
    PECGratingReducedML(PeriodicCurve* curve, double angle, double wavenumber, std::vector<unsigned> numElementsVector, std::shared_ptr<std::vector<EmpiricalInterpolation>> empiricalGreen, std::shared_ptr<std::vector<EmpiricalInterpolation>> empiricalRhs, std::shared_ptr<BEM::Matrix> bemPodMatrix, std::shared_ptr<std::vector<BEM::Matrix>> bemMatrices, std::shared_ptr<std::vector<BEM::ColVector>> bemRhs, unsigned greenTerms = 50U);
    ~PECGratingReducedML(void);
    void buildDiscrete(void) override {};
    int solve(void) override {return 0;};
    double matrixError(void) override;
    const BEM::Matrix &interpolatedMatrix(void) override;
private:
    void interpolateGreen(void);
    std::vector<unsigned> _numElementsVector;
    std::shared_ptr<std::vector<EmpiricalInterpolation>> _empiricalGreen = nullptr;
    std::shared_ptr<std::vector<EmpiricalInterpolation>> _empiricalRhs = nullptr;
    std::shared_ptr<BEM::Matrix> _bemPodMatrix = nullptr;
    std::shared_ptr<std::vector<BEM::Matrix>> _bemMatrices = nullptr;
    std::shared_ptr<std::vector<BEM::ColVector>> _bemRhs = nullptr;
    std::vector<std::unique_ptr<DiscreteSpaceOnCurve_1D>> _spaces;
    std::vector<std::unique_ptr<Operator>> _operators;
};



class UQProblemFactoryReducedML : public UQProblemFactoryReduced
{
public:
    UQProblemFactoryReducedML(double period, double angle, double wavenumber, double pertSize, int dimDecay, std::vector<unsigned> numElementsVector, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients, unsigned greenTerms = 50U);
    ~UQProblemFactoryReducedML();
    // std::unique_ptr<Problem> generateNewProblem(unsigned counter) const override {return nullptr;};
    void trainEmpirical([[maybe_unused]] unsigned samples, [[maybe_unused]] double tolerance = 1E-6) override {};
    void trainEmpirical([[maybe_unused]] unsigned samples, [[maybe_unused]] unsigned bases) override {};
    void trainEmpirical(const std::vector<unsigned> &samples, unsigned bases);
    std::unique_ptr<PECGrating> generateProblemFromCurve(PeriodicCurve *curve) override;
protected:
    void saveEIData(void) override;
    void sampleMC(const std::vector<unsigned> &samples);
    void initializeSVD(void) override;
    PECGratingML generateNewTrainingProblem(unsigned level, unsigned counter);
    std::unique_ptr<Problem> generateNewProblem(unsigned counter) const override;
    std::vector<unsigned> _numElementsVector;
    std::shared_ptr<std::vector<EmpiricalInterpolation>> _empiricalGreen = nullptr;
    std::shared_ptr<std::vector<EmpiricalInterpolation>> _empiricalRhs = nullptr;
    std::shared_ptr<std::vector<EmpiricalInterpolation>> _empiricalSolution = nullptr;
    std::shared_ptr<std::vector<BEM::Matrix>> _bemPodBasis = nullptr;
    std::shared_ptr<std::vector<BEM::Matrix>> _bemMatrices = nullptr;
    std::shared_ptr<std::vector<BEM::Matrix>> _bemRedMatrices = nullptr;
    std::shared_ptr<std::vector<BEM::ColVector>> _bemRhs = nullptr;
    std::shared_ptr<std::vector<BEM::ColVector>> _bemRedRhs = nullptr;
    std::shared_ptr<BEM::Matrix> _bemPodMatrix = nullptr;
private:
    std::vector<std::unique_ptr<DiscreteSpaceOnCurve_1D>> _spaces;
};



#endif
