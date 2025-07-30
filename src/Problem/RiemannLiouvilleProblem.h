#ifndef RIE_LIO_PROB
#define RIE_LIO_PROB

#include <ParametricProblem.h>
#include <MyTypes.h>
#include <memory>
#include <vector>
#include <unordered_map> 
#include <Curve.h>
#include <Mesh.h>
#include <functional>
#include <RiemannLiouville.h>

//////////////////////////////////////////////////////////////
// Classes for Riemann Liouville problem derivative and RB. //
//////////////////////////////////////////////////////////////


/**
 * @brief: Class for managing structures regarding Greedy RB construction.
 * TODO: Move to its own file.
 */
class GreedyHelper
{
public:
    GreedyHelper(BEM::Matrix prodMatrix, std::vector<BEM::Matrix> matrixVector, std::vector<BEM::ColVector> rhsVector, const std::function<double (std::vector<double>)> &infSupEst);
    void loadBasis(BEM::ColVector newBasis);

    BEM::ColVector solveAtPoint(std::vector<double> point);
    BEM::ColVector solveAtPointReduced(std::vector<double> point);
    double errorAtPoint(std::vector<double> point);
private:
    BEM::Matrix reducedMatrixAtPoint(std::vector<double> point);
    BEM::ColVector reducedRhsAtPoint(std::vector<double> point);
    BEM::Matrix _pMat;
    Eigen::ColPivHouseholderQR<BEM::Matrix> _qr;
    std::unique_ptr<BEM::Matrix> _basis=nullptr;
    std::vector<BEM::Matrix> _hifiMatrices;
    std::vector<BEM::Matrix> _reducedMatrices;
    std::vector<BEM::ColVector> _hifiRhs;
    std::vector<BEM::ColVector> _reducedRhs;
    const std::function<double (std::vector<double>)> &_infSupEst;
    // for residual estimation
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, BEM::Complex>>>> _MatMatProd;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, BEM::Complex>>> _MatRhsProd;
    std::unordered_map<int, std::unordered_map<int, BEM::Complex>> _RhsRhsProd;
};

// Forward declarations
class ExplicitScalarFunction_1D;
class DiscreteSpaceOnCurve_1D;
class SobolevSlobodeckii;
class Operator;
class Curve2D;
class Mesh1D;
class MeshCurve1D;
class DiscreteSpaceMesh; 
class ExplicitScalarFunction_2D;
class LeftFracDerivative;
class DiscreteFunctionMesh;
using VectorFun_1D = std::vector<ExplicitScalarFunction_1D>;
using VectorFun_2D = std::vector<ExplicitScalarFunction_2D>;
class ScalarFunctionBase_1D;
class ScalarFunctionBase_2D;

/**
 * @brief: Class for Riemman liouville problem.
 */
class RiemannLiouvilleProblem : public ProblemMesh
{
public:
    RiemannLiouvilleProblem(int order, DiscreteSpaceMesh &space);
    template <typename T, typename S, typename R>
    RiemannLiouvilleProblem(int order, DiscreteSpaceMesh &space, T qFun, S dFun, R fFun) :
        RiemannLiouvilleProblem(order, space)
    {
        _rhsFun.reset(new R(fFun));
        _RL.reset(new RiemannLiouvilleMesh(_space, RiemannLiouvilleMesh::Side::LEFT, _order, dFun, qFun));
    }

    RiemannLiouvilleProblem(int order, DiscreteSpaceMesh &space, std::shared_ptr<ScalarFunctionBase_1D> &qFun, std::shared_ptr<ScalarFunctionBase_1D> &dFun, std::shared_ptr<ScalarFunctionBase_2D> &fFun) :
        RiemannLiouvilleProblem(order, space)
    {
        _rhsFun = fFun;
        _RL.reset(new RiemannLiouvilleMesh(_space, RiemannLiouvilleMesh::Side::LEFT, _order, dFun, qFun));
    }

    ~RiemannLiouvilleProblem(void);
    void buildDiscrete(void) override;
    const DiscreteSpaceMesh &getSpace(void) override {return _space;}
private:
    int _order;
    const DiscreteSpaceMesh &_space;
    const Mesh1D &_mesh;
    std::shared_ptr<ScalarFunctionBase_2D> _rhsFun = nullptr;
    std::unique_ptr<Operator> _RL = nullptr;
};

/**
 * @brief: Factory for Riemann liouville problems and greedy training.
 */
class RiemannLiouvilleMeshFactory : public ProblemFactory
{
public:
    RiemannLiouvilleMeshFactory(unsigned numberOfElements, int order, std::vector<ExplicitScalarFunction_1D> qVector, std::vector<ExplicitScalarFunction_1D> dVector, std::vector<ExplicitScalarFunction_2D> fVector);
    std::unique_ptr<ProblemMesh> generateNewProblem(std::vector<double> parameters);
    void trainGreedy(std::vector<BEM::Interval1D> limits, int points, double tolerance, const std::function<double (std::vector<double>)> &infSupEst, const std::vector<std::vector<double>> *testingPoints=nullptr);
    BEM::ColVector greedySolve(std::vector<double> parameters);
    BEM::Matrix getSobMat(void);
    ~RiemannLiouvilleMeshFactory(void);
private:
    unsigned _numberOfElements = 0u;    
    int _order = 0;
    StraightCurve _curve;
    MeshCurve1D _mesh;
    size_t _dimensionQ;
    size_t _dimensionD;
    size_t _dimensionF;
    std::vector<std::shared_ptr<ExplicitScalarFunction_1D>> _qVector;
    std::vector<std::shared_ptr<ExplicitScalarFunction_1D>> _dVector;
    std::vector<std::shared_ptr<ExplicitScalarFunction_2D>> _fVector;
    std::unique_ptr<GreedyHelper> _greedy = nullptr;
    std::unique_ptr<LeftFracDerivative> _LF;
    std::unique_ptr<DiscreteSpaceMesh> _space = nullptr;
};

#endif 
