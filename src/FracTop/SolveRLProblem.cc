
#include <SolveRLProblem.h>
#include <Curve.h>
#include <Utilities.h>
#include <memory>
#include <RegularP1Mesh.h>
#include <ScalarValuedFunction.h>
#include <ParametrizedFunction.h>
#include <RiemannLiouvilleProblem.h>
#include <RegularP1Mesh.h>

void solveRLProblem(int order, const std::string &typeD, const std::string &typeQ, const BEM::CVector &dVec, const BEM::CVector &qVec, const BEM::CVector &rhsVec, int ms)
{
    std::shared_ptr<ScalarFunctionBase_1D> dFun;
    std::shared_ptr<ScalarFunctionBase_1D> qFun;
    if (typeD == "poly") {
        dFun.reset(new PolynomialFunction_1D(dVec));
    } else if (typeD == "trig") {
        dFun.reset(new TrigonometricFunction_1D(1, dVec));
    } else if (typeD == "pw") {
        dFun.reset(new PwConstantFunction_1D(0, 1, dVec));
    }

    if (typeQ == "poly") {
        qFun.reset(new PolynomialFunction_1D(qVec));
    } else if (typeQ == "trig") {
        qFun.reset(new TrigonometricFunction_1D(1, qVec));
    } else if (typeQ == "pw") {
        qFun.reset(new PwConstantFunction_1D(0, 1, qVec));
    }
    PolynomialFunction_1D rhsFun(rhsVec);
    ExplicitScalarFunction_2D rhs2D([&rhsFun](double t, double s [[maybe_unused]]){return rhsFun(t);});
        
    TrigonometricCurve curve(1, 0, std::vector<double>{0}, std::vector<double>{0});
    MeshCurve1D mesh(ms, curve);
    RegularP1_0Mesh_1D space(mesh);
    RiemannLiouvilleProblem RLP(order, space, qFun, dFun, rhs2D);

    RLP.buildDiscrete();
    RLP.solve();

    auto sol = RLP.getSolutionVec();
    auto solution = space.generateFunction(BEM::toVector(sol));
    auto &solutionDerL = solution->derivative(order);
    BEM::plotFunction("SOL", *solution, mesh);
    BEM::plotFunction("SOLDER", solutionDerL, mesh);
    BEM::plotFunction("DIF", *dFun);
    BEM::plotFunction("REC", *qFun);
    BEM::plotFunction("RHS", rhsFun);
}