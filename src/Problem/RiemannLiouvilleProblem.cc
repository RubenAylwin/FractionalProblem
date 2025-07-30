#include <RiemannLiouvilleProblem.h>
#include <SobolevSlobodeckii.h>
#include <L2Operator.h>
#include <DiscreteSpace.h>
#include <ScalarValuedFunction.h>
#include <Curve.h>
#include <GreenQP.h>
#include <Operator.h>
#include <OneDimensionalIntegration.h>
#include <TwoDimensionalIntegration.h>
#include <Msg.h>
#include <Utilities.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <RegularP1Mesh.h>
#include <FractionalDerivative.h>

useMessages("RIE_PROB");

/**
 * @brief: Constructor. Takes in all matrices for different parameters.
 */
GreedyHelper::GreedyHelper(BEM::Matrix prodMatrix, std::vector<BEM::Matrix> matrixVector, std::vector<BEM::ColVector> rhsVector, const std::function<double (std::vector<double>)> &infSupEst) :
    _pMat{prodMatrix},
    _basis(nullptr),
    _hifiMatrices{matrixVector},
    _hifiRhs{rhsVector},
    _infSupEst{infSupEst}
{
    _qr = _pMat.colPivHouseholderQr();
    for (size_t i = 0; i < _hifiRhs.size(); ++i) {
        for (size_t j = 0; j < _hifiRhs.size(); ++j) {
            _RhsRhsProd[i][j] = _hifiRhs[i].dot(_qr.solve(_hifiRhs[j]));
        }
    }
}

/**
 * @brief: Loading a new basis and updates the matrices and rhs.
 */
void GreedyHelper::loadBasis(BEM::ColVector newBasis)
{
    msg(5) << "start GreedyHelper::loadBasis" << endMsg;
    if (not _basis) {
        _basis.reset(new BEM::Matrix(newBasis.size(), 1));
        _basis->col(0) = newBasis;
    } else {
        _basis->conservativeResize(_basis->rows(), _basis->cols() + 1);
        _basis->col(_basis->cols() - 1) = newBasis;
    }
    _reducedMatrices.clear();
    for (const auto &matrix : _hifiMatrices) {
        _reducedMatrices.push_back(BEM::product(BEM::product(_basis->adjoint(), matrix), *_basis));
    }
    _reducedRhs.clear();
    for (const auto &rhs : _hifiRhs) {
        _reducedRhs.push_back(BEM::product(_basis->adjoint(), rhs));
    }

    for (size_t i = 0; i < _hifiMatrices.size(); ++i) {
        for (size_t j = 0; j < _hifiMatrices.size(); ++j) {
            for (int k = 0; k < _basis->cols(); ++k) {
                for (int l = 0; l < _basis->cols(); ++l) {
                    _MatMatProd[i][j][k][l] = (_hifiMatrices[i]*_basis->col(k)).dot(_qr.solve(_hifiMatrices[j]*_basis->col(l)));
                }
            }
        }
    }

    for (size_t i = 0; i < _hifiMatrices.size(); ++i) {
        for (int j = 0; j < _basis->cols(); ++j) {
            for (size_t k = 0; k < _hifiRhs.size(); ++k) {
                _MatRhsProd[i][j][k] = (_hifiMatrices[i]*_basis->col(j)).dot(_qr.solve(_hifiRhs[k]));
            }
        }
    }

    msg(5) << "end GreedyHelper::loadBasis" << endMsg;
}

/**
 * @brief: Given a point, return the reduced matrix for the parameters.
 */
BEM::Matrix GreedyHelper::reducedMatrixAtPoint(std::vector<double> point)
{
    msg(5) << "start GreedyHelper::reducedMatrixAtPoint" << endMsg;
    return BEM::linearCombination<BEM::Matrix>(_reducedMatrices, std::vector<double>{point.begin(), point.begin() + _reducedMatrices.size()});
    msg(5) << "end GreedyHelper::reducedMatrixAtPoint" << endMsg;
}

/**
 * @brief: Given a point, return the reduced RHS for the parameters.
 */
BEM::ColVector GreedyHelper::reducedRhsAtPoint(std::vector<double> point)
{
    return BEM::linearCombination<BEM::ColVector>(_reducedRhs, std::vector<double>{point.begin() + _reducedMatrices.size(), point.end()});
}

/**
 * @brief: Solve the problem at a point. Returns the solution expressed in the high fidelity basis.
 */
BEM::ColVector GreedyHelper::solveAtPoint(std::vector<double> point)
{
    msg(5) << "start GreedyHelper::solveAtPoint" << endMsg;
    if (_reducedMatrices.size() == 0) {
        BEM::ColVector result = BEM::ColVector::Zero(_hifiRhs[0].size());
        return result;
    }
    BEM::ColVector result = *_basis*solveAtPointReduced(point);
    msg(5) << "end GreedyHelper::solveAtPoint" << endMsg;
    return result;
}

/**
 * @brief: Solve the problem at a point. Returns expression in the reduced basis.
 */
BEM::ColVector GreedyHelper::solveAtPointReduced(std::vector<double> point)
{
    msg(5) << "start GreedyHelper::solveAtPointReduced" << endMsg;
    BEM::Matrix redMat = reducedMatrixAtPoint(point);
    msg(5) << "GreedyHelper::solveAtPointReduced got redMat" << endMsg;
    BEM::ColVector redRHS = reducedRhsAtPoint(point);
    msg(5) << "GreedyHelper::solveAtPointReduced got redRHS" << endMsg;
    BEM::ColVector redSol = redMat.colPivHouseholderQr().solve(redRHS);
    msg(5) << "GreedyHelper::solveAtPointReduced got Sol" << endMsg;
    msg(5) << "end GreedyHelper::solveAtPointReduced" << endMsg;
    return redSol;
}

/**
 * @brief: ecomputes the (estimated) error at a point.
 */
double GreedyHelper::errorAtPoint(std::vector<double> point)
{
    // Different error components.
    BEM::Complex errorMatMat = 0.0;
    BEM::Complex errorMatRhs = 0.0;
    BEM::Complex errorRhsRhs = 0.0;
    
    // Only compute error from basis if there is one.
    // Important for the begining where no elements have been chosen yet.
    if (_basis) {
        BEM::ColVector redSol = solveAtPointReduced(point);
        for (size_t i = 0; i < _hifiMatrices.size(); ++i) {
            for (size_t j = 0; j < _hifiMatrices.size(); ++j) {
                for (int k = 0; k < _basis->cols(); ++k) {
                    for (int l = 0; l < _basis->cols(); ++l) {
                        errorMatMat += redSol[k]*redSol[l]*point[i]*point[j]*_MatMatProd[i][j][k][l];
                    }
                }
            }
        }

        for (size_t i = 0; i < _hifiMatrices.size(); ++i) {
            for (int j = 0; j < _basis->cols(); ++j) {
                for (size_t k = 0; k < _hifiRhs.size(); ++k) {
                    errorMatRhs -= 2.*redSol[j]*point[i]*point[_hifiMatrices.size() + k]*_MatRhsProd[i][j][k];
                }
            }
        }
    }
    for (size_t i = 0; i < _hifiRhs.size(); ++i) {
        for (size_t j = 0; j < _hifiRhs.size(); ++j) {
            errorRhsRhs += point[_hifiMatrices.size() + i]*point[_hifiMatrices.size() + j]*_RhsRhsProd[i][j];
        }
    }

    return std::sqrt(std::abs(errorMatMat + errorMatRhs + errorRhsRhs))/(_infSupEst(point));
}

/**
 * @brief: Problem constructor.
 * @desc: order is to be given as a number between 50 and 100. The real fractional order is order/100.
 * This is necessary because fractional derivatives are saved in an int-indexed dictionary with the same format.
 */
RiemannLiouvilleProblem::RiemannLiouvilleProblem(int order, DiscreteSpaceMesh &space) :
    _order{order},
    _space{space},
    _mesh{space.getMesh()}
{
}

/**
 * @brief: Explicit destructor.
 */
RiemannLiouvilleProblem::~RiemannLiouvilleProblem()
{
}

/**
 * @brief: Build the discrete problem.
 */
void RiemannLiouvilleProblem::buildDiscrete(void)
{
    msg(5) << "RiemannLiouvilleQ::buildDiscrete Start" << endMsg;
    assert(_RL);
    _RL->assembleMassMatrix();
    msg(5) << "RiemannLiouvilleQ::buildDiscrete done with RL" << endMsg;
    _matrix.reset(new BEM::Matrix(_RL->getMatrix()));
    _rhs.reset(new BEM::ColVector(_space.testAgainstBasis(*_rhsFun)));
    msg(5) << "RiemannLiouvilleQ::buildDiscrete End" << endMsg;
}

/**
 * @brief: Factory constructor. Same considerations as before for the order.
 */
RiemannLiouvilleMeshFactory::RiemannLiouvilleMeshFactory(unsigned numberOfElements, int order, std::vector<ExplicitScalarFunction_1D> qVector, std::vector<ExplicitScalarFunction_1D> dVector, std::vector<ExplicitScalarFunction_2D> fVector) :
    _numberOfElements(numberOfElements),
    _order(order),
    _curve(1.),
    _mesh(_numberOfElements, _curve),
    _dimensionQ{qVector.size()},
    _dimensionD{dVector.size()},
    _dimensionF{fVector.size()},
    _qVector(qVector.size(), nullptr),
    _dVector(dVector.size(), nullptr),
    _fVector(fVector.size(), nullptr),
    _space{new RegularP1_0Mesh_1D(_mesh)}
{
    msg(5) << "start RiemannLiouvilleMeshFactory::Constructor" << endMsg;
    for (size_t i = 0; i < _dimensionQ; ++i) {
        _qVector[i].reset(new ExplicitScalarFunction_1D(qVector[i]));
    }
    for (size_t i = 0; i < _dimensionD; ++i) {
        _dVector[i].reset(new ExplicitScalarFunction_1D(dVector[i]));
    }
    for (size_t i = 0; i < _dimensionF; ++i) {
        _fVector[i].reset(new ExplicitScalarFunction_2D(fVector[i]));
    }
    msg(5) << "end RiemannLiouvilleMeshFactory::Constructor" << endMsg;
}

/**
 * @brief: Generate a new problem for the given parameters.
 */
std::unique_ptr<ProblemMesh> RiemannLiouvilleMeshFactory::generateNewProblem(std::vector<double> parameters)
{
    msg(5) << "start RiemannLiouvilleMeshFactory::generateNewProblem" << endMsg;
    auto dQ = _dimensionQ;
    auto dD = _dimensionD;
    auto dF = _dimensionF;
    assert(parameters.size() == dQ + dD + dF);
    
    ExplicitScalarFunction_1D q([&, parameters, dQ](double t) -> BEM::Complex {
        BEM::Complex result = 0.0;
        for (size_t i = 0; i < dQ; ++i) {
            result += _qVector[i]->operator()(t)*parameters[i];
        }
        return result;
    });

    ExplicitScalarFunction_1D d([&, parameters, dD, dQ](double t) -> BEM::Complex {
        BEM::Complex result = 0.0;
        for (size_t i = dQ; i < dQ + dD; ++i) {
            result += _dVector[i - dQ]->operator()(t)*parameters[i];
        }
        return result;
    });

    ExplicitScalarFunction_2D f([&, parameters, dF, dQ, dD](double t, double s) -> BEM::Complex {
        BEM::Complex result = 0.0;
        for (size_t i = dQ + dD; i < dQ + dD + dF; ++i) {
            result += _fVector[i - dQ - dD]->operator()(t, s)*parameters[i];
        }
        return result;
    });
    msg(5) << "end RiemannLiouvilleMeshFactory::generateNewProblem" << endMsg;    
    return std::unique_ptr<ProblemMesh>(new RiemannLiouvilleProblem(_order, *_space, q, d, f));
}

/**
 * @brief: Get a RB.
 */
void RiemannLiouvilleMeshFactory::trainGreedy(std::vector<BEM::Interval1D> limits, int points, double tolerance, const std::function<double (std::vector<double>)> &infSupEst, const std::vector<std::vector<double>> *testingPoints)
{
    msg(5) << "start RiemannLiouvilleMeshFactory::trainGreedy" << endMsg;
    assert(limits.size() == _dimensionQ + _dimensionD + _dimensionF);
    assert(tolerance > 1E-9);

    //Points
    std::vector<std::vector<double>> quadPoints;
    // std::vector<double> firstPoint{};
    {
        std::vector<std::vector<double>> toTensPoints(_dimensionQ + _dimensionD + _dimensionF, std::vector<double>());
        GaussLegendre_1D integ(points);

        for (size_t i = 0; i < _dimensionQ + _dimensionD + _dimensionF; ++i) {
            // If a dimension is given same lower and upper limit, then we dont tensorize there.
            // Bit of a hacky solution to being able to give functions that will have no real
            // parameter on the parametrized problem.
            if (limits[i].first == limits[i].second) {
                toTensPoints[i] = std::vector<double>{limits[i].first};
            } else {
                toTensPoints[i] = integ.points(limits[i].first, limits[i].second);
            }
            // firstPoint.push_back(i < _dimensionQ ? 0.0 : 1.0);
        }
        quadPoints = BEM::tensorize(toTensPoints);
    }
    
    // HiFi Rhs
    std::vector<BEM::ColVector> completeRhs{};

    // HiFi matrices for affine decomposition
    std::vector<BEM::Matrix> completeMatrices{};

    // auto nominalProblem = generateNewProblem(firstPoint);
    // nominalProblem->buildDiscrete();
    // nominalProblem->solve();
    // auto nominalSolution = BEM::ColVector(nominalProblem->getSolutionVec());

    int counter = 0;
    // Save matrices for affine decomposition
    for (size_t i = 0; i < _dimensionQ + _dimensionD; ++i) {
        std::vector<double> auxPoint(_dimensionQ + _dimensionD + _dimensionF, 0.0);
        auxPoint[i]=1.0;
        auto prob = generateNewProblem(auxPoint);
        prob->buildDiscrete();
        completeMatrices.push_back(prob->getMatrix());
    }

    // Save RHS for affine decomposition
    for (size_t i = 0; i < _dimensionF; ++i) {
        std::vector<double> auxPoint(_dimensionQ + _dimensionD + _dimensionF, 0.0);
        auxPoint[_dimensionQ + _dimensionD + i] = 1.0;
        auto prob = generateNewProblem(auxPoint);
        prob->buildDiscrete();
        completeRhs.push_back(prob->getRHS());
    }

    // Get matrix for duality product. Here we are using leftDer*leftDer.
    _LF.reset(new LeftFracDerivative(*_space, _order));
    _LF->assembleMassMatrix();
    auto lfMat = _LF->getMatrix();

    // Create a greedy helper with the generated data
    _greedy.reset(new GreedyHelper(lfMat, completeMatrices, completeRhs, infSupEst));
    GreedyHelper &greedy = *_greedy;

    
    // maxError in samples (initialized so that we enter the while loop)
    double maxError = tolerance + 1;

    while (maxError > tolerance) {
        maxError = 0.0;
        std::vector<double> toAdd;
        for (const auto &point : quadPoints) {
            double error = greedy.errorAtPoint(point);
            if (error > maxError) {
                maxError = error;
                toAdd = point;
            }
        }
        BEM::ColVector addedSolution = BEM::linearCombination(completeMatrices, std::vector<double>(toAdd.begin(), toAdd.begin() + _dimensionQ + _dimensionD)).colPivHouseholderQr().solve(BEM::linearCombination(completeRhs, std::vector<double>(toAdd.begin() + _dimensionQ + _dimensionD, toAdd.end())));
        if (maxError > tolerance) {
            auto realSolution = _space->generateFunction(addedSolution);
            auto &realSolutionDer = realSolution->derivative(_order);
            levelIf(1) {
                BEM::ColVector residual = BEM::linearCombination(completeMatrices, std::vector<double>(toAdd.begin(), toAdd.begin() + _dimensionQ + _dimensionD))*greedy.solveAtPoint(toAdd) - BEM::linearCombination(completeRhs, std::vector<double>(toAdd.begin() + _dimensionQ + _dimensionD, toAdd.end()));
                double realError = std::sqrt(std::abs(residual.dot(lfMat.colPivHouseholderQr().solve(residual))));
                msg(1) << "Max Error on " << counter << ": " << maxError << " - " << realError << " at point " << endMsg;
                if(counter){
                    auto greedySolVec = greedy.solveAtPoint(toAdd);
                    auto greedySolution = _space->generateFunction(greedySolVec);
                    auto &greedySolutionDer = greedySolution->derivative(_order);
                    realSolutionDer -= greedySolutionDer;
                    msg(1) << "real Error L2Norm: " << realSolutionDer.L2Norm() << endMsg;
                    BEM::plotFunction("solAtPoint"+std::to_string(counter), *greedySolution, _mesh);
                } else {
                    msg(1) << "real Error L2Norm: " << realSolutionDer.L2Norm() << endMsg;
                }
                BEM::plotFunction("basis"+std::to_string(counter), *(_space->generateFunction(addedSolution)), _mesh);
            }
            greedy.loadBasis(addedSolution/addedSolution.norm());
            counter++;
        } else {
            levelIf(1) {
                msg(1) << "Done with max error on " << counter << ": " << maxError << endMsg;
                auto greedySolution = _space->generateFunction(greedy.solveAtPoint(toAdd));
                auto &greedySolutionDer = greedySolution->derivative(_order);
                auto realSolution = _space->generateFunction(addedSolution);
                auto &realSolutionDer = realSolution->derivative(_order);
                realSolutionDer -= greedySolutionDer;
                msg(1) << "real Error: " << realSolutionDer.L2Norm() << endMsg;
                BEM::plotFunction("solAtPoint"+std::to_string(counter), *greedySolution, _mesh);
            }
        }
        if (testingPoints) {
            double maxTestError = 0.0;
            BEM::ColVector maxTestSolution;
            BEM::ColVector maxGreedySolution;
            for (auto &point : *testingPoints) {
                BEM::ColVector testSolution =  BEM::linearCombination(completeMatrices, std::vector<double>(point.begin(), point.begin() + _dimensionQ + _dimensionD)).colPivHouseholderQr().solve(BEM::linearCombination(completeRhs, std::vector<double>(point.begin() + _dimensionQ + _dimensionD, point.end())));
                auto greedySolution = greedy.solveAtPoint(point);
                double error = (testSolution - greedySolution).norm()/std::sqrt(_numberOfElements);
                if (error > maxTestError) {
                    maxTestError = error;
                    maxTestSolution = testSolution;
                    maxGreedySolution = greedySolution;
                }
            }
            std::cout << "Max test Error := " << maxTestError << std::endl;
        }
    }
}

/**
 * @brief: Explicit destructor.
 */
RiemannLiouvilleMeshFactory::~RiemannLiouvilleMeshFactory(void)
{
}

/**
 * @brief: Solve with the greedy problem.
 */
BEM::ColVector RiemannLiouvilleMeshFactory::greedySolve(std::vector<double> parameters)
{
    return _greedy->solveAtPoint(parameters);
}

/**
 * @brief: Get the matrix for norm-product.
 */
BEM::Matrix RiemannLiouvilleMeshFactory::getSobMat(void)
{
    return _LF->getMatrix();
}
