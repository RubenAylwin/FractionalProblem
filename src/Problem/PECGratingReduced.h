#ifndef PEC_GRATING_REDUCED
#define PEC_GRATING_REDUCED

#include <ParametricProblem.h>
#include <PECGrating.h>
#include <MyTypes.h>
#include <memory>
#include <vector>
#include <random>

///////////////////////////////////////////////////////////
// Classes refering to RB versions of PECGrating problem //
///////////////////////////////////////////////////////////

//Forward declaration
class EmpiricalInterpolation;

/**
 * @brief: Reduced problem. Requires an empirical interpolator.
 */
class PECGratingReduced : public PECGrating
{
public:
    PECGratingReduced(PeriodicCurve* curve, double angle, double wavenumber, unsigned numElements, const BEM::EmpiricalInterpolationPtr &empiricalGreen, const BEM::EmpiricalInterpolationPtr &empiricalRhs,  const std::shared_ptr<std::vector<BEM::Matrix>> &greenMatrices,  const std::shared_ptr<std::vector<BEM::ColVector>> &rhsVectors , const std::shared_ptr<BEM::Matrix> &bemPodBasis, unsigned greenTerms = 50U);
    ~PECGratingReduced();
    void buildDiscrete(void) override;
    int solve(void) override;
    double matrixError(void) override;
    const BEM::Matrix &interpolatedMatrix(void) override;
    const BEM::ColVector &getRHS(void) override;
    double rhsError(void);
    double solutionError(void);
private:
    BEM::EmpiricalInterpolationPtr _empiricalGreen = nullptr;
    BEM::EmpiricalInterpolationPtr _empiricalRhs = nullptr;
    std::shared_ptr<std::vector<BEM::Matrix>> _greenMatrices = nullptr;
    std::shared_ptr<std::vector<BEM::ColVector>> _rhsVectors = nullptr;
    std::shared_ptr<BEM::Matrix> _bemPodBasis = nullptr;
};

/**
 * @brief: Reduced Factory. First, it generates training problems, then it builds a RB using SVN and generates reduced problems.
 */
class UQProblemFactoryReduced : public UQProblemFactory
{
public:
    UQProblemFactoryReduced(double period, double angle, double wavenumber, double pertSize, int dimDecay, unsigned numElements, std::vector<double> sineCoefficients, std::vector<double> cosineCoefficients, unsigned greenTerms = 50U);
    virtual void trainEmpirical(unsigned samples, double tolerance = 1E-6);
    virtual void trainEmpirical(unsigned samples, unsigned bases);
    ~UQProblemFactoryReduced();
    virtual std::unique_ptr<PECGrating> generateProblemFromCurve(PeriodicCurve *curve);
protected:
    std::unique_ptr<Problem> generateNewProblem(unsigned counter) const override;
    // std::unique_ptr<Problem> generateNewTestingProblem(unsigned counter) const override;
    void sampleMC(unsigned samples);
    virtual void initializeSVD(void);
    virtual void saveEIData(void);
    std::shared_ptr<std::vector<BEM::Matrix>> _greenMatrices = nullptr;
    std::shared_ptr<std::vector<BEM::ColVector>> _rhsVectors = nullptr;
    bool _trained = false;
    std::shared_ptr<BEM::Matrix> _bemPodBasis = nullptr;
    
private:
    std::shared_ptr<EmpiricalInterpolation> _empiricalGreen = nullptr;
    std::shared_ptr<EmpiricalInterpolation> _empiricalRhs = nullptr;
    std::shared_ptr<EmpiricalInterpolation> _empiricalSolution = nullptr;
};


#endif 
