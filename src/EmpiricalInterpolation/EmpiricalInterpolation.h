#ifndef EMPIRICAL_INTERPOLATION
#define EMPIRICAL_INTERPOLATION
#include <MyTypes.h>
#include <vector>
#include <list>

///////////////////////////////////////////////
// Class helper for Empirical Interpolation. //
///////////////////////////////////////////////

class EmpiricalInterpolation {
public:
    EmpiricalInterpolation(unsigned dimension);
    EmpiricalInterpolation(std::string name, unsigned dimension);
    ~EmpiricalInterpolation() = default;
    void loadTrainingData(std::vector<BEM::Complex> &&data);
    void loadTrainingData(BEM::ColVector &&data);
    void buildInterpolator(double tolerance);
    void buildInterpolator(unsigned bases);
    void computeSVD(void);
    const BEM::ColVector &getSingValues(void) const {return _singVals;}
    const std::vector<int> &getInterpolationIndices(void) const;
    BEM::ColVector interpolate(const std::vector<BEM::Complex> &indices) const;
    BEM::ColVector interpolate(const BEM::ColVector &indices) const;
    BEM::ColVector getSolution(const std::vector<BEM::Complex> &indices) const;
    BEM::ColVector getSolution(const BEM::ColVector &indices) const;
    const BEM::Matrix &podBasis(void) const {return _podBasis;}
protected:
    std::string _name = "";
    const unsigned _dimension = 0;
    std::vector<BEM::ColVector> _trainingData;
    std::vector<int> _indices;
    BEM::Matrix _podBasis;
    BEM::Matrix _projMatrix;
    BEM::ColVector _singVals;
    std::shared_ptr<BEM::Matrix> _U = nullptr;
};

#endif
