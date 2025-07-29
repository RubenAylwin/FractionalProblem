#ifndef GENERAL_PROBLEM
#define GENERAL_PROBLEM

#include <memory>
#include <vector>
#include <MyTypes.h>
#include <chrono>
class GreenHQP2D;
class QPProblem;
class DiscreteSpaceOnCurve_1D;
class DiscreteSpaceMesh;

////////////////////////////////////////
// Base class for parametric problems //
////////////////////////////////////////

class Problem {
public:
    Problem(void) = default;
    virtual ~Problem(void) = default;
    virtual void buildDiscrete(void) = 0;
    virtual const DiscreteSpaceOnCurve_1D &getSpace(void) = 0;
    virtual const BEM::ColVector &getSolutionVec(void) { return *_solution; };
    virtual const BEM::ColVector &getRHS(void) { return *_rhs; };
    virtual const BEM::Matrix &getMatrix(void) const { return *_matrix; };
    virtual BEM::Complex getMatrixIndex(int row, int col) const;
    virtual BEM::Complex getRHSIndex(int row) const;
    virtual int solve(void);
    virtual QPProblem *toQP(void) const { return nullptr; };
protected:
    std::unique_ptr<BEM::Matrix> _matrix = nullptr;
    std::unique_ptr<BEM::ColVector> _solution = nullptr;
    std::unique_ptr<BEM::ColVector> _rhs = nullptr;
};

class ProblemMesh {
public:
    ProblemMesh(void) = default;
    virtual ~ProblemMesh(void) = default;
    virtual void buildDiscrete(void) = 0;
    virtual const DiscreteSpaceMesh &getSpace(void) = 0;
    virtual const BEM::ColVector &getSolutionVec(void) { return *_solution; };
    virtual const BEM::ColVector &getRHS(void) { return *_rhs; };
    virtual const BEM::Matrix &getMatrix(void) const { return *_matrix; };
    virtual BEM::Complex getMatrixIndex(int row, int col) const;
    virtual BEM::Complex getRHSIndex(int row) const;
    virtual int solve(void);
protected:
    std::unique_ptr<BEM::Matrix> _matrix = nullptr;
    std::unique_ptr<BEM::ColVector> _solution = nullptr;
    std::unique_ptr<BEM::ColVector> _rhs = nullptr;
};


class ProblemFactory {
public:
    ProblemFactory(void) = default;
    virtual ~ProblemFactory(void) = default;
};

class RandomProblemFactory : public ProblemFactory {
public:
    RandomProblemFactory(void) = default;
    virtual void reset(void) { _counter = 0U; };
    virtual std::unique_ptr<Problem> generateNewProblem(void);
    virtual std::unique_ptr<Problem> generateNewTestingProblem(void);
protected:
    virtual std::unique_ptr<Problem> generateNewProblem(unsigned counter) const = 0;
    virtual std::unique_ptr<Problem> generateNewTestingProblem(unsigned counter) const = 0;
private:
    unsigned _counter = 0U;
    unsigned _testingCounter = 0U;

};

///////////////////////////////////
// Specialization to QP problems //
///////////////////////////////////

class QPProblem : public Problem {
public:
    QPProblem(double period, double angle, double wavenumber);
    virtual ~QPProblem(void);
    void setGreenWindowTerms(unsigned terms);
    QPProblem *toQP(void) const override { return const_cast<QPProblem *>(this); };
protected:
    double _period = 0;
    double _angle = 0;
    double _wavenumber = 0;
    std::unique_ptr<GreenHQP2D> _green = nullptr;
};

class QPProblemFactory : public ProblemFactory {
public:
    QPProblemFactory(unsigned numElements);
    virtual ~QPProblemFactory() {};
    unsigned getSize(void) {return _numElements;}
protected:
    unsigned _numElements = 0;
};

class QPProblemFactoryRB : public QPProblemFactory {
public:
    QPProblemFactoryRB(unsigned numElements);
    virtual ~QPProblemFactoryRB() {};
    virtual std::unique_ptr<Problem> generateNewProblem2(unsigned counter, bool empirical = false) const = 0;
};


#endif
