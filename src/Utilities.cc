#include <Utilities.h>
#include <MyTypes.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ScalarValuedFunction.h>
#include <Curve.h>
#include <boost/math/special_functions/hankel.hpp>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include <Eigen/Eigenvalues>
#include <EmpiricalInterpolation.h>
#include <Msg.h>
#include <DiscreteSpace.h>
#include <unordered_set>
#include <set>

useMessages("UTILS_MSG");

BEM::Complex BEM::cexp(BEM::Complex t)
{
    return std::exp(t);
}

BEM::Complex BEM::cexp(double t)
{
    return BEM::Complex(std::cos(t), std::sin(t));
}

std::unique_ptr<ScalarFunctionBase_1D> BEM::generatePlaneWave(const double angle, const double wavenumber, const Curve2D &curve)
{
    auto planeWave = [=, &curve](const double t)->BEM::Complex {
        return cexp(wavenumber*(std::sin(angle)*curve.at(t).getX() - std::cos(angle)*curve.at(t).getY()))*curve.jacobian(t);
    };
    return std::unique_ptr<ScalarFunctionBase_1D>(new ExplicitScalarFunction_1D(planeWave));
}

std::unique_ptr<ScalarFunctionBase_1D> BEM::generatePlaneWave(const BEM::Complex xFactor, BEM::Complex yFactor, const Curve2D &curve)
{

    auto planeWave = [=, &curve](const double t)->BEM::Complex {
        return cexp((xFactor*curve.at(t).getX() + yFactor*curve.at(t).getY()))*curve.jacobian(t);
    };
    return std::unique_ptr<ScalarFunctionBase_1D>(new ExplicitScalarFunction_1D(planeWave));
}


std::unique_ptr<ScalarFunctionBase_2D> BEM::generatePlaneWave(const double angle, const double wavenumber)
{
    auto planeWave = [=](const double x, const double y)->BEM::Complex {
        return cexp(wavenumber*(std::sin(angle)*x - std::cos(angle)*y));
    };
    return std::unique_ptr<ScalarFunctionBase_2D>(new ExplicitScalarFunction_2D(planeWave));
}


int BEM::maxVectorIndex(const BEM::ColVector &vec)
{
    int index = 0;
    double max = 0.0;
    for (unsigned row = 0; row < vec.rows(); ++row) {
        if (max < std::abs(vec[row])) {
            max = std::abs(vec[row]);
            index = row;
        }
    }
    
    return index;
}

double BEM::inftyError(const BEM::Matrix &original, const BEM::Matrix &approximation, bool relative)
{
    double max = 0.0;
    for (unsigned row = 0; row < original.rows(); ++row) {
        for (unsigned col = 0; col < original.cols(); ++col) {
            double value = relative ? std::abs(original(row, col) - approximation(row, col))/std::abs(original(row, col)) : std::abs(original(row, col) - approximation(row, col));
            if (max < value) {
                max = value;
            }
        }
    }
    
    return max;
}


double BEM::pError(const BEM::Matrix &original, const BEM::Matrix &approximation, bool relative)
{
    return relative ? (original - approximation).norm()/original.norm() : (original - approximation).norm();
}


template<int N>
struct A {
    constexpr A() : arr() {
        arr[0] = 1;
        for (auto i = 1; i != N; ++i)
            arr[i] = arr[i-1]*(-(2.*i - 1)*(2.*i - 1)/(i*8.));
    }
    double arr[N];
};

[[maybe_unused]] constexpr static auto AGR = A<200>();

[[maybe_unused]] static BEM::Complex expansion(int order, double argument)
{
    BEM::Complex result = 0;
    BEM::Complex term = 1;
    for (int i = 0; i < order; ++i) {
        result += term*AGR.arr[i];
        term *= BEM::I/(argument);
        if (std::abs(term/result) <= 1E-6) {
            break;
        }
    }
    return result;
}

BEM::Complex BEM::hankel0_1(double z)
{
    return boost::math::cyl_hankel_1(0, z);
}

BEM::Complex BEM::hankel0_1_orig(double z)
{
    return boost::math::cyl_hankel_1(0, z);   
}

BEM::ColVector BEM::conjugate(const ColVector &vec) {
    BEM::ColVector result(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i) {
        result[i] = std::conj(vec[i]);
    }
    return result;
}

BEM::ColVector BEM::toColVector(const std::vector<BEM::Complex> &vector)
{
    BEM::ColVector vec(vector.size());
    for (unsigned i = 0; i < vector.size(); ++i) {
        vec[i] = vector[i];
    }

    return vec;
}

BEM::ColVector BEM::toColVector(const std::vector<double> &vector)
{
    BEM::ColVector vec(vector.size());
    for (unsigned i = 0; i < vector.size(); ++i) {
        vec[i] = vector[i];
    }

    return vec;
}

std::vector<BEM::Complex> BEM::toVector(const BEM::ColVector &vector)
{
    std::vector<BEM::Complex> vec(vector.size());
    for (unsigned i = 0; i < vector.size(); ++i) {
        vec[i] = vector[i];
    }

    return vec;
}


BEM::ColVector BEM::stretch(const BEM::Matrix &matrix)
{
    BEM::ColVector result(matrix.rows()*matrix.cols());
    int i = 0;
    for (unsigned col = 0; col < matrix.cols(); ++col) {
        for (unsigned row = 0; row < matrix.rows(); ++row) {
            result[i] = matrix(row, col);
            ++i;
        }
    }
    assert(i == matrix.rows()*matrix.cols());
    return result;
}

BEM::Matrix BEM::compress(const BEM::ColVector &vector)
{
    double N = 0;
    assert(std::modf(std::sqrt(vector.size()), &N) == 0.0);
    assert(N*N == vector.size());
    BEM::Matrix matrix((int)N, (int)N);
    int i = 0;
    for (unsigned col = 0; col < matrix.cols(); ++col) {
        for (unsigned row = 0; row < matrix.rows(); ++row) {
            matrix(row, col) = vector[i];
            ++i;
        }
    }
    assert(i == matrix.rows()*matrix.cols());
    return matrix;
}

BEM::Interval1D BEM::intersect(const Interval1D &sup1, const Interval1D &sup2)
{
    double a = std::max(sup1.first, sup2.first);
    double b = std::min(sup1.second, sup2.second);
    return Interval1D(a, b);
}

template<>
bool BEM::getEnv(std::string envVar)
{
    auto var = std::getenv(envVar.c_str());

    if (not var) {
        return false;
    }
    
    if (std::string(var) == "true") {
        return true;
    }

    return false;
}

template<>
int BEM::getEnv(std::string envVar)
{
    auto var = std::getenv(envVar.c_str());

    if (not var) {
        return -1;
    }
    
    return std::stoi(var);
}

void BEM::squareDiagonal(BEM::Matrix &diag)
{
    assert(diag.rows() == diag.cols());
    for (unsigned i = 0; i < diag.rows(); ++i) {
        if (std::abs(diag(i, i)) > 0) {
            diag(i, i) *= diag(i, i);
        }
    }
}

void BEM::sqrtDiagonal(BEM::Matrix &diag)
{
    assert(diag.rows() == diag.cols());
    for (unsigned i = 0; i < diag.rows(); ++i) {
        if (std::abs(diag(i, i)) > 0) {
            diag(i, i) = std::sqrt(diag(i, i));
        }
    }
}

void BEM::invertDiagonal(BEM::Matrix &diag)
{
    assert(diag.rows() == diag.cols());
    for (unsigned i = 0; i < diag.rows(); ++i) {
        if (std::abs(diag(i, i)) > 0) {
            diag(i, i) = 1.0/diag(i, i);
        }
    }
}

BEM::Matrix BEM::transpose(const BEM::Matrix &matrix)
{

    BEM::Matrix tMatrix(matrix.cols(), matrix.rows());
    tbb::parallel_for( tbb::blocked_range2d<int>(0, matrix.cols(), 0, matrix.rows()),
                       [&]( const tbb::blocked_range2d<int> &r ) {
                           for(int j=r.rows().begin(), j_end=r.rows().end(); j<j_end; ++j){
                               for(int k=r.cols().begin(), k_end=r.cols().end(); k<k_end; ++k){
                                   tMatrix(j, k) = matrix(k, j);
                               }
                           }
                       });

    return tMatrix;
    
}

BEM::Matrix BEM::product(const BEM::Matrix &matrix1, const BEM::Matrix &matrix2)
{
    assert(matrix1.cols() == matrix2.rows());
    BEM::Matrix pMatrix(matrix1.rows(), matrix2.cols());
    tbb::parallel_for( tbb::blocked_range2d<int>(0, matrix1.rows(), 0, matrix2.cols()),
                       [&]( const tbb::blocked_range2d<int> &r ) {
                           for(int j=r.rows().begin(), j_end=r.rows().end(); j<j_end; ++j){
                               for(int k=r.cols().begin(), k_end=r.cols().end(); k<k_end; ++k){
                                   pMatrix(j, k) = matrix1.row(j).conjugate().dot(matrix2.col(k));
                               }
                           }
                       });
    return pMatrix;
}

double BEM::distance(const Interval1D &sup1, const Interval1D &sup2)
{
    auto intersection = intersect(sup1, sup2);

    if (intersection.first <= intersection.second) {
        return 0.;
    }

    return intersection.first - intersection.second;
        
}

BEM::TimePoint BEM::now(void)
{
    return std::chrono::system_clock::now();
}

int BEM::timeDifference(const BEM::TimePoint &time1, const BEM::TimePoint &time2)
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
}

void BEM::getBasesPerLevel(std::vector<EmpiricalInterpolation> &empiricalInterpolators, unsigned bases)
{
    msg(5) << "(start)BEM::getBasesPerLevel" << endMsg;
    msg(6) << "Called with " << bases << " bases" << endMsg;
    std::vector<unsigned> basesPerLevel(empiricalInterpolators.size(), 0U);
    for (unsigned b = 0; b < bases; ++b) {
        int largest = -1;
        double value = 0;
        for (unsigned level = 0; level < empiricalInterpolators.size(); ++level) {
            double largestAtLevel = basesPerLevel[level] < empiricalInterpolators[level].getSingValues().size() ? std::pow(std::abs(empiricalInterpolators[level].getSingValues()[basesPerLevel[level]]), 2) : 0;
            if (value < largestAtLevel) {
                value = largestAtLevel;
                largest = level;
            }
        }
        if (largest >=0) {
            ++basesPerLevel[largest];
        }
    }
    msg(5) << "(mid)BEM::getBasesPerLevel" << endMsg;    
    for (unsigned level = 0; level < empiricalInterpolators.size(); ++level) {
        if (basesPerLevel[level] == 0) {
            ++basesPerLevel[level];
        }
        msg(6) << "(mid)BEM::getBasesPerLevel level " << level << " has " << basesPerLevel[level] << " bases" << endMsg;
        empiricalInterpolators[level].buildInterpolator(basesPerLevel[level]);
        msg(6) << "built level " << level << endMsg;
    }
    msg(5) << "(end)BEM::getBasesPerLevel" << endMsg;
}

void BEM::getBasesPerLevel(std::vector<EmpiricalInterpolation> &empiricalInterpolators, double tolerance)
{
    double norm{0.0};
    for (auto &emp : empiricalInterpolators) {
        emp.computeSVD();
        for (const auto &val : emp.getSingValues()) {
            norm += std::abs(val)*std::abs(val);
        }
    }
    double energy{0};
    std::vector<unsigned> bases(empiricalInterpolators.size(), 0U);
    while (energy/norm < 1. - std::pow(tolerance, 2)) {
        unsigned largest = 0;
        double value = 0;
        for (unsigned level = 0; level < empiricalInterpolators.size(); ++level) {
            double largestAtLevel = bases[level] < empiricalInterpolators[level].getSingValues().size() ? std::pow(std::abs(empiricalInterpolators[level].getSingValues()[bases[level]]), 2) : 0;
            if (value < largestAtLevel) {
                value = largestAtLevel;
                largest = level;
            }
        }
        ++bases[largest];
        energy += value;
        
    }

    for (unsigned level = 0; level < empiricalInterpolators.size(); ++level) {
        empiricalInterpolators[level].buildInterpolator(bases[level]);
    }

}
/**
 * @brief: Given a partition, find the partition index it lies in. If it lies on a partition node, the given direction
 * indicates in which interval to place the point (LEFT->limit at point from the left).
 */
unsigned BEM::findIndex(const std::vector<double> &partition, double point, BasisFunction_1D::Direction direction)

{
    unsigned lowerIndex = 0;
    unsigned higherIndex = partition.size() - 1;
    while (lowerIndex + 1 < higherIndex) {
        unsigned newIndex = std::floor(0.5*(lowerIndex + higherIndex));
        if (partition[newIndex] < point) {
            lowerIndex = newIndex;
        } else if (partition[newIndex] > point) {
            higherIndex = newIndex;
        } else if (partition[newIndex] == point) {
            switch(direction) {
            case (BasisFunction_1D::Direction::LEFT):
                return newIndex > 0 ? newIndex - 1 : newIndex;
                break;
            case (BasisFunction_1D::Direction::RIGHT):
                return newIndex < partition.size() - 1 ? newIndex : newIndex - 1;
            }
        }
    }
    return lowerIndex;
}

std::vector<BEM::Complex> BEM::basisVector(int i, int size)
{
    std::vector<BEM::Complex> result(size, 0.0);
    result[i] = 1.0;
    return result;
}

double BEM::L2Norm(const ScalarFunctionBase_1D &function)
{
    int numPoints = 2000;
    double result = 0;
    for(double j = 1; j < numPoints; ++j) {
        result += std::real(function(j/numPoints)*function(j/numPoints))/(numPoints-1);
    }
    return std::sqrt(result);
}

BEM::Support1DL BEM::join(const BEM::Support1DL &sup1, const BEM::Support1DL &sup2)
{
    std::unordered_set<double> includedSet;
    BEM::Support1DL result;
    for (const auto &val : sup1) {
        result.emplace_back(val);
        assert(includedSet.insert(val.first).second);
    }
    for (const auto &val : sup2) {
        if(includedSet.insert(val.first).second) {
            result.emplace_back(val);
        }
    }
    return result;
}

BEM::Support1DL BEM::remove(const BEM::Support1DL &sup1, const BEM::Support1DL &sup2)
{
    std::unordered_set<double> remove;
    BEM::Support1DL result;
    for (const auto &val : sup1) {
        assert(remove.insert(val.first).second);
    }
    for (const auto &val : sup2) {
        if(remove.find(val.first) == remove.end()) {
            result.emplace_back(val);
        }
    }
    return result;
}


void BEM::plotFunction(std::string fileName, const ScalarFunctionBase_1D &function)
{
    int numPoints = 200;
    BEM::DVector absVal, realVal, imagVal, position;
    std::ofstream myFile;
    myFile.open(fileName+".txt");
    for(double j = 0; j < numPoints; ++j) {
        position.push_back(j/numPoints);
        absVal.push_back(std::abs(function(j/numPoints)));
        realVal.push_back(function(j/numPoints).real());
        imagVal.push_back(function(j/numPoints).imag());
        myFile << position.back() << " " << function(j/numPoints).real() << " " << function(j/numPoints).imag() << std::endl;
        msg(5) << position.back() << " - " << function(j/numPoints) << endMsg;
    }
    myFile.close();
}

void BEM::plotFunction(std::string fileName, const DiscreteFunctionMesh &function, const Mesh1D &mesh)
{
    int numPoints = 20;
    BEM::DVector absVal, realVal, imagVal, position;
    std::ofstream myFile;
    myFile.open(fileName+".txt");
    for (unsigned index = 0; index < mesh.numElements(); ++index) {
        auto element = mesh.getElement(index);
        for (double j = 0; j < numPoints; ++j) {
            double innerPoint = (1.*j)/numPoints;
            position.push_back(element(innerPoint).getX());
            absVal.push_back(std::abs(function.evaluate(index, j/numPoints)));
            realVal.push_back(function.evaluate(index, j/numPoints).real());
            imagVal.push_back(function.evaluate(index, j/numPoints).imag());
            myFile << position.back() << " " << function.evaluate(index, j/numPoints).real() << " " << function.evaluate(index, j/numPoints).imag() << std::endl;
            msg(5) << position.back() << " - " << function.evaluate(index, j/numPoints) << endMsg;
        }
    }
    myFile.close();
}

std::vector<std::vector<double>> BEM::tensorize(const std::vector<std::vector<double>> points)
{
    std::vector<size_t> position(points.size(), 0);
    int totalPoints = 1;
    for (const auto &vec : points) {
        totalPoints *= vec.size();
    }
    std::vector<std::vector<double>> result;
    for (int p = 0; p < totalPoints; ++p) {
        std::vector<double> currentPoint(points.size(), 0);
        for (size_t i = 0; i < points.size(); ++i) {
            currentPoint[i] = points[i][position[i]];
        }
        result.push_back(currentPoint);        
        for (size_t i = 0; i < points.size(); ++i) {
            position[i] += 1;
            if (position[i] < points[i].size()) {
                break;
            }
            position[i] = 0;
        }
    }
    return result;
}

std::vector<BEM::Matrix> BEM::reduceMatrices(const std::vector<BEM::Matrix> &fullMatrices, const BEM::Matrix &reducedBasis)
{
    std::vector<BEM::Matrix> reducedMatrices{};
    for (const auto &matrix : fullMatrices) {
        reducedMatrices.push_back(BEM::product(BEM::product(reducedBasis.adjoint(), matrix), reducedBasis));
    }
    return reducedMatrices;
}

