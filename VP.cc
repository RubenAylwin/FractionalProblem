#include <Curve.h>
#include <DiscreteSpace.h>
#include <ScalarValuedFunction.h>
#include <ParametrizedFunction.h>
#include <MyTypes.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <Utilities.h>
#include <chrono>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <unordered_set>
#include <progressbar.hpp>
#include <RegularP1Mesh.h>
#include <RiemannLiouville.h>
#include <RiemannLiouvilleProblem.h>
#include <Msg.h>
useMessages("MAIN");
namespace po = boost::program_options;


template <typename T>
T readOption(po::variables_map &vm, std::string &&option, T &&defVal) {
    T value;
    try {
        value = vm[option].as<T>();
    } catch (std::exception &e) {
        std::cout << "Caught exception when reading --" << option <<": ";
        std::cout << e.what() << std::endl;
        value = defVal;
    }
    return value;
}

void solveRLProblem(int order, std::string typeD, std::string typeQ, BEM::CVector dVec, BEM::CVector qVec, BEM::CVector rhsVec, int ms)
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
    ExplicitScalarFunction_2D rhs2D([&rhsFun](double t, double s){return rhsFun(t);});
        
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


void evaluateRBforRL(int order, std::string typeD, std::string typeQ, BEM::CVector dVec, BEM::CVector qVec, BEM::CVector dVarVec, BEM::CVector qVarVec, BEM::CVector rhsVec, int ms)
{
    std::vector<BEM::Interval1D> limits{};
    
    VectorFun_1D qFunVec{};
    for (size_t i = 0; i < qVec.size(); ++i) {
        BEM::CVector auxQVec(qVec.size(), 0);
        auxQVec[i] = 1.;
        limits.push_back(BEM::Interval1D(std::real(qVec[i] - qVarVec[i]), std::real(qVec[i] + qVarVec[i])));
        if (typeQ == "poly") {
            PolynomialFunction_1D qAux(auxQVec);
            qFunVec.push_back(ExplicitScalarFunction_1D([qAux](double t){return qAux(t);}));
        } else if (typeQ == "trig") {
            TrigonometricFunction_1D qAux(1, auxQVec);
            qFunVec.push_back(ExplicitScalarFunction_1D([qAux](double t){return qAux(t);}));
        } else if (typeQ == "pw") {
            PwConstantFunction_1D qAux(0, 1, auxQVec);
            qFunVec.push_back(ExplicitScalarFunction_1D([qAux](double t){return qAux(t);}));
        }
        BEM::plotFunction("Q"+std::to_string(i), qFunVec.back());
    }

    VectorFun_1D dFunVec{};
    for (size_t i = 0; i < dVec.size(); ++i) {
        BEM::CVector auxDVec(dVec.size(), 0.0);
        limits.push_back(BEM::Interval1D(std::real(dVec[i] - dVarVec[i]), std::real(dVec[i] + dVarVec[i])));
        auxDVec[i] = 1.;
        if (typeD == "poly") {
            PolynomialFunction_1D dAux(auxDVec);
            dFunVec.push_back(ExplicitScalarFunction_1D([dAux](double t){return dAux(t);}));
        } else if (typeD == "trig") {
            TrigonometricFunction_1D dAux(1, auxDVec);
            dFunVec.push_back(ExplicitScalarFunction_1D([dAux](double t){return dAux(t);}));
        } else if (typeD == "pw") {
            PwConstantFunction_1D dAux(0, 1, auxDVec);
            dFunVec.push_back(ExplicitScalarFunction_1D([dAux](double t){return dAux(t);}));
        }
        BEM::plotFunction("D"+std::to_string(i), dFunVec.back());
    }

    PolynomialFunction_1D rhsFun(rhsVec);
    ExplicitScalarFunction_2D rhs2D([&rhsFun](double t, double s){return rhsFun(t);});
    limits.push_back(BEM::Interval1D(1., 1.));

    
    auto est = [order](std::vector<double> point) -> double {
        // double min = *std::min_element(point.begin() + 2, point.begin() + 4 + 2);
        // double max = *std::max_element(point.begin() + 2, point.begin() + 4 + 2);
        return std::abs(std::cos(M_PI*order/100.));
        // return 0.5*((max + min)*std::abs(std::cos(M_PI*fOrder)) - (max - min));
    };
    RiemannLiouvilleMeshFactory factory(ms, order, qFunVec, dFunVec, VectorFun_2D{rhs2D});
    factory.trainGreedy(limits, 100, 1e-6, est);
}

int main(int argc, char* argv[]) {
    po::options_description desc("Options");
    desc.add_options()("Problem", po::value<int>(), "Specify problem to solve")
        ("help", "Show help menu")
        ("rhs", po::value<BEM::CVector>()->multitoken(), "Taylor expansion for right hand side.")
        ("dt", po::value<std::string>(), "Expansion type for diffusion coefficient.")
        ("qt", po::value<std::string>(), "Expansion type for reaction coefficient.")
        ("d", po::value<BEM::CVector>()->multitoken(), "Expansion for diffusion coefficient.")
        ("q", po::value<BEM::CVector>()->multitoken(), "Expansion for reaction coefficient.")
        ("dv", po::value<BEM::CVector>()->multitoken(), "Variation for diffusion coefficient.")
        ("qv", po::value<BEM::CVector>()->multitoken(), "Variation for reaction coefficient.")
        ("fo", po::value<int>(), "Fractional order.")
        ("mesh_size", po::value<int>(), "Num. mesh elements.");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, ~po::command_line_style::allow_short), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    if (not vm.count("dt") or not vm.count("qt")){
        std::cerr << "Missing data type (--dt and/or --qt)." << std::endl;
        std::cerr << "Quitting without solving." << std::endl;
        return -1;
    }

    auto dVec = readOption(vm, std::string("d"), BEM::CVector{1.});
    auto qVec = readOption(vm, std::string("q"), BEM::CVector{0.});
    auto rhsVec = readOption(vm, std::string("rhs"), BEM::CVector{1.});
    auto typeD = vm["dt"].as<std::string>();
    auto typeQ = vm["qt"].as<std::string>();
    int ms = readOption(vm, std::string("mesh_size"), 100);
    int order = readOption(vm, std::string("fo"), 70);
            
    if (vm["Problem"].as<int>() == 1){
        std::cout << "Solving Riemann Liouville Problem with specified parameters." << std::endl;    
        solveRLProblem(order, typeD, typeQ, dVec, qVec, rhsVec, ms);
        return 0;
    } else if (vm["Problem"].as<int>() == 2){

        auto dVarVec = readOption(vm, std::string("dv"), BEM::CVector(dVec.size(), 0.0));
        auto qVarVec = readOption(vm, std::string("qv"), BEM::CVector(qVec.size(), 0.0));

        if (dVarVec.size() != dVec.size()){
            std::cerr << "Data --d and --dv must have the same length" << std::endl;
            std::cerr << "Quitting without solving." << std::endl;
            return -1;
        }

        if (qVarVec.size() != qVec.size()){
            std::cerr << "Data --q and --qv must have the same length" << std::endl;
            std::cerr << "Quitting without solving." << std::endl;
            return -1;
        }
        std::cout << "Constructing and evaluating a Reduced Basis for the Riemann Liouville Problem." << std::endl;
        evaluateRBforRL(order, typeD, typeQ, dVec, qVec, dVarVec, qVarVec, rhsVec, ms);
    } else {
        std::cerr << "Invalid problem type." << std::endl;
        std::cerr << "Quitting without solving." << std::endl;
        return -1;
    }
    return 0;
}


// int main2(int argc, char* argv[]) {
//     po::options_description desc("Options");
//     desc.add_options()("Problem", po::value<int>(), "Specify problem to solve")
//         ("help", "Show help menu")
//         ("P", po::value<double>(), "Period")
//         ("WL", po::value<double>(), "Wavelength")
//         ("WN", po::value<double>(), "Wavenumber")
//         ("ANG", po::value<double>(), "Angle")
//         ("DOF", po::value<int>(), "NumOfElements")
//         ("OK", po::value<int>(), "NumOfElements for Overkill")
//         ("SIN", po::value<std::vector<double>>()->multitoken(), "SinCoef")
//         ("COS", po::value<std::vector<double>>()->multitoken(), "CosCoef")
//         ("DCY", po::value<double>(), "Decay")
//         ("PS", po::value<double>(), "Perturbation size")
//         ("TG", po::value<int>(), "TrainingSamples")
//         ("TT", po::value<int>(), "TestingSamples")
//         ("RL", po::value<int>(), "Realizations")
//         ("RBG", po::value<int>(), "RBGreen");
    
//     po::variables_map vm;
//     po::store(po::parse_command_line(argc, argv, desc), vm);
//     po::notify(vm);
//     if (vm.count("help")) {
//         std::cout << desc << std::endl;
//         return 1;
//     }
    
//     if (vm.count("P")) {
//         if (vm["Problem"].as<int>() == 1){
//             runAndPlotSolutionForHolo(vm["P"].as<double>(), vm.count("WL") ? 2.*M_PI/vm["WL"].as<double>() : vm["WN"].as<double>(),  vm["ANG"].as<double>(), vm["DOF"].as<int>(), vm["SIN"].as<std::vector<double>>(), vm["COS"].as<std::vector<double>>());
//             return 0;
//         }

//         if (vm["Problem"].as<int>() == 2) {
//             SLEmpiricalTest(vm["P"].as<double>(), vm.count("WL") ? 2.*M_PI/vm["WL"].as<double>() : vm["WN"].as<double>(),  vm["ANG"].as<double>(), vm["SIN"].as<std::vector<double>>(), vm["COS"].as<std::vector<double>>(), vm["PS"].as<double>(), vm["DCY"].as<double>(), vm["DOF"].as<int>(), vm["OK"].as<int>(), vm["TG"].as<int>(), vm["TT"].as<int>(), vm["RBG"].as<int>(), vm["RL"].as<int>());
//         }

//         if (vm["Problem"].as<int>() == 3) {
//             Reliability(vm["P"].as<double>(), vm.count("WL") ? 2.*M_PI/vm["WL"].as<double>() : vm["WN"].as<double>(),  vm["ANG"].as<double>(), vm["SIN"].as<std::vector<double>>(), vm["COS"].as<std::vector<double>>(), vm["PS"].as<double>(), vm["DCY"].as<double>(), vm["DOF"].as<int>(), vm["TG"].as<int>(), vm["RBG"].as<int>());
//         }
        
//     }
//     return 0;
// }
