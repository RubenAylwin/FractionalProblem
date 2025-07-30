#ifndef UTILITIES
#define UTILITIES
    
#include <complex>
#include <memory>
#include <MyTypes.h>
#include <string>
#include <vector>
#include <DiscreteSpace.h>
#include <DiscreteSpaceMesh.h>
#include <Mesh.h>
#include <chrono>

#define startTimer(timerName) auto start_##timerName = std::chrono::system_clock::now();
#define stopTimer(timerName, message) auto stop_##timerName = std::chrono::system_clock::now(); \
    msg(0) << "Time taken for " << message << ": " << std::chrono::duration_cast<std::chrono::microseconds>(stop_##timerName - start_##timerName).count() << "[ms]" << endMsg;

class ScalarFunctionBase_1D;
class Curve2D;
class EmpiricalInterpolation;

namespace BEM {
    
    std::complex<double> cexp(BEM::Complex t);
    std::complex<double> cexp(double t);
    std::unique_ptr<ScalarFunctionBase_1D> generatePlaneWave(const double angle, const double wavenumber, const Curve2D &curve);
    std::unique_ptr<ScalarFunctionBase_2D> generatePlaneWave(const double angle, const double wavenumber);
    std::unique_ptr<ScalarFunctionBase_1D> generatePlaneWave(const BEM::Complex xFactor, BEM::Complex yFactor, const Curve2D &curve);
    int maxVectorIndex(const ColVector &vec);
    double inftyError(const Matrix &original, const Matrix &approximation,  bool relative = false);
    double pError(const Matrix &original, const Matrix &approximation,  bool relative = false);
    Complex hankel0_1(double z);
    Complex hankel0_1_orig(double z);
    ColVector conjugate(const ColVector &vec);
    ColVector toColVector(const std::vector<BEM::Complex> &vector);
    ColVector toColVector(const std::vector<double> &vector);
    std::vector<BEM::Complex> toVector(const ColVector &vector);
    ColVector stretch(const BEM::Matrix &matrix);
    Matrix compress(const BEM::ColVector &vector);
    Interval1D intersect(const Interval1D &sup1, const Interval1D &sup2);
    template<typename T>
    T getEnv(std::string envVar);
    void invertDiagonal(BEM::Matrix &diag);
    void squareDiagonal(BEM::Matrix &diag);
    void sqrtDiagonal(BEM::Matrix &diag);
    BEM::Matrix transpose(const BEM::Matrix &matrix);
    BEM::Matrix product(const BEM::Matrix &matrix1, const BEM::Matrix &matrix2);
    double distance(const Interval1D &sup1, const Interval1D &sup2);
    BEM::TimePoint now(void);
    int timeDifference(const BEM::TimePoint &time1, const BEM::TimePoint &time2);
    void getBasesPerLevel(std::vector<EmpiricalInterpolation> &empiricalInterpolators, unsigned bases);
    void getBasesPerLevel(std::vector<EmpiricalInterpolation> &empiricalInterpolators, double tolerance);
    unsigned findIndex(const std::vector<double> &partition, double point, BasisFunction_1D::Direction direction);
    std::vector<BEM::Complex> basisVector(int i, int size);
    double L2Norm(const ScalarFunctionBase_1D &function);
    Support1DL join(const Support1DL &sup1, const Support1DL &sup2);
    Support1DL remove(const BEM::Support1DL &sup1, const BEM::Support1DL &sup2);
    void plotFunction(std::string fileName, const ScalarFunctionBase_1D &function);
    void plotFunction(std::string fileName, const DiscreteFunctionMesh &function, const Mesh1D &mesh);
    std::vector<std::vector<double>> tensorize(const std::vector<std::vector<double>> points);
    std::vector<BEM::Matrix> reduceMatrices(const std::vector<BEM::Matrix> &fullMatrices, const BEM::Matrix &reducedBasis);
    template <typename T, typename S>
    T affineCombination(std::vector<T> vec, S point)
    {
        assert(vec.size() == point.size() + 1);
        T result = vec[0];
        for (size_t index = 0; index < point.size(); ++index) {
            result += point[index]*vec[index+1];
        }
        return result;
    }
    template <typename T, typename S>
    T linearCombination(std::vector<T> vec, S point)
    {
        assert(vec.size() == point.size());
        T result = vec[0]*point[0];
        for (size_t index = 1; index < point.size(); ++index) {
            result += point[index]*vec[index];
        }
        return result;
    }
    template <typename T>
    void printVector(const std::vector<T> &vec) {
        std::cout << "[" << vec[0];
        for (int i = 1; i < vec.size(); ++i) {
            std::cout << ", " << vec[i];
        }
        std::cout << "]\n";
    }
    template <typename T>
    T vectorNorm(const std::vector<T> &vec) {
        T result = std::abs(vec[0]*vec[0]);
        for (int i = 1; i < vec.size(); ++i) {
            result += std::abs(vec[i]*vec[i]);
        }
        return std::sqrt(result);
    }

}

#endif
