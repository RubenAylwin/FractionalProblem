#ifndef OPERATOR
#define OPERATOR

#include <MyTypes.h>
#include <memory>

///////////////////////////////////
// Base class for all operators. //
///////////////////////////////////

class Operator {
public:
    Operator(const unsigned testSize, const unsigned trialSize);
    virtual ~Operator() = default;
    virtual BEM::Complex indexedDuality(const unsigned i, const unsigned j) = 0;
    virtual const BEM::Matrix &getMatrix(void);
    virtual void assembleMassMatrix(void);
protected:
    std::unique_ptr<BEM::Matrix> _massMatrix;
    const unsigned _testSize;
    const unsigned _trialSize;
};

#endif
