#include <ParametricProblem.h>
#include <Curve.h>
#include <Utilities.h>
#include <WeaklySingular.h>
#include <Halton.h>
#include <GreenQP.h>
#include <DiscreteSpace.h>
#include <iostream>
#include <fstream>
#include <MyTypes.h>


/**
 * @brief: egenerate a new problem. Updates the counter.
 */
std::unique_ptr<Problem> RandomProblemFactory::generateNewProblem(void)
{
    _counter++;
    return generateNewProblem(_counter);
}

/**
 * @brief: Generates new testing problem. Updates different counter.
 */
std::unique_ptr<Problem> RandomProblemFactory::generateNewTestingProblem(void)
{
    _testingCounter++;
    return generateNewTestingProblem(_testingCounter);
}

/**
 * @brief: Get matrix info.
 */
BEM::Complex Problem::getMatrixIndex(int row, int col) const
{
    assert(_matrix);
    return (*_matrix)(row, col);
}

/**
 * @brief: Return the index of the RHS.
 */
BEM::Complex Problem::getRHSIndex(int row) const
{
    assert(_rhs);
    return (*_rhs)[row];
}

/**
 * @brief: Solve the discrete problem.
 */
int Problem::solve(void)
{
    auto first = BEM::now();
    _solution.reset(new BEM::ColVector(_matrix->colPivHouseholderQr().solve(*_rhs)));
    auto second = BEM::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(second - first).count();
}

/**
 * @brief: Get Matrix info.
 */
BEM::Complex ProblemMesh::getMatrixIndex(int row, int col) const
{
    assert(_matrix);
    return (*_matrix)(row, col);
}

/**
 * @brief: Return the index of the RHS.
 */
BEM::Complex ProblemMesh::getRHSIndex(int row) const
{
    assert(_rhs);
    return (*_rhs)[row];
}

/**
 * @brief: Solve the discrete problem.
 */
int ProblemMesh::solve(void)
{
    auto first = BEM::now();
    _solution.reset(new BEM::ColVector(_matrix->colPivHouseholderQr().solve(*_rhs)));
    auto second = BEM::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(second - first).count();
}

/**
 * @brief: Constructor.
 */
QPProblem::QPProblem(double period, double angle, double wavenumber) :
    _period{period},
    _angle{angle},
    _wavenumber{wavenumber},
    _green(new GreenHQP2D(period, angle, wavenumber))
{
}

/**
 * @brief: Destructor.
 */
QPProblem::~QPProblem(void)
{
}

/**
 * @brief: Access to Green's function computation
 */
void QPProblem::setGreenWindowTerms(unsigned terms)
{
    _green->setWindowTerms(terms);
}

/**
 * @brief: Constructor.
 */
QPProblemFactory::QPProblemFactory(unsigned numElements) :
    _numElements{numElements}
{
}

/**
 * @brief: Constructor.
 */
QPProblemFactoryRB::QPProblemFactoryRB(unsigned numElements) :
    QPProblemFactory{numElements}
{
}

