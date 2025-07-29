#include <EmpiricalInterpolation.h>
#include <Utilities.h>
#include <numeric>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <chrono>
#include <Msg.h>
#include <progressbar.hpp>
#include <algorithm>
#include "randomized_svd.h"

useMessages("EMP_INT");


/**
 * @brief: Constructor.
 * @desc: If given a name, something will be different about how the sing vals are saved.
 */
EmpiricalInterpolation::EmpiricalInterpolation(std::string name, unsigned dimension) :
    _name(name),
    _dimension(dimension),
    _trainingData(),
    _indices(),
    _singVals(),
    _U{nullptr}
{
    _trainingData.reserve(200);
}

/**
 * @brief: Constructor without name.
 */
EmpiricalInterpolation::EmpiricalInterpolation(unsigned dimension) :
    _dimension(dimension),
    _trainingData(),
    _indices(),
    _singVals(),
    _U{nullptr}
{
}

/**
 * @brief: Add a new column to the training data.
 */
void EmpiricalInterpolation::loadTrainingData(BEM::ColVector &&data)
{
    //Correctness check
    assert(data.size() == _dimension);
    _trainingData.push_back(std::move(data));
}

/**
 * @brief: Overload with std::vector<>
 */
void EmpiricalInterpolation::loadTrainingData(std::vector<BEM::Complex> &&data)
{
    loadTrainingData(BEM::toColVector(data));
}

/**
 * @brief: Compute the SVD of the data matrix.
 */
void EmpiricalInterpolation::computeSVD(void)
{
    if (_trainingData.size() == 0) {
        return;
    }
    
    BEM::Matrix dataMatrix(_dimension, _trainingData.size());
    for (unsigned col = 0; col < dataMatrix.cols(); ++col) {
        for (unsigned row = 0; row < dataMatrix.rows(); ++row) {
            dataMatrix(row, col) = _trainingData[col](row);
        }
    }
    _trainingData.clear();
    auto start = std::chrono::system_clock::now();
    
    // If the NO_SVD env var is active, then we do not do the SVD and take the matrix as is
    if (BEM::getEnv<bool>("NO_SVD")) {
        msg(0) << "Working without SVD" << endMsg;
        _U.reset(new BEM::Matrix(dataMatrix));
        return;
    }

    // If the RAND_SVD env var is active, then we use a random version of the SVD
    if (BEM::getEnv<bool>("RAND_SVD")) {
        msg(0) << "Working with random_svd" << endMsg;
        RandomizedSvd rsvd(dataMatrix, 400, 100, 1); // Hardcoded values for now.
        _U.reset(new BEM::Matrix(rsvd.matrixU()));
        _singVals = rsvd.singularValues();
    } else {
        Eigen::BDCSVD<BEM::Matrix, Eigen::ComputeThinU> SVD(dataMatrix);
        _singVals = SVD.singularValues();
        _U.reset(new BEM::Matrix(SVD.matrixU()));
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    msg(0) << "SVD Time: " << elapsed.count() << "[s]" << std::endl;
    msg(5) << "The singular values are: " << endMsg;
    msg(5) << _singVals << endMsg;
}

/**
 * @brief: build the empirical interpolation with a given number of bases.
 */
void EmpiricalInterpolation::buildInterpolator(unsigned bases)
{
    msg(5) << "(start)EmpiricalInterpolation::buildInterpolator" << endMsg;
    msg(5) << "Called to generate " << bases << " bases" << endMsg;
    // If SVD is not computed, compute it.
    if (not _U) {
        assert(_trainingData.size() > 0 and bases > 0);
        computeSVD();
    }

    _podBasis = BEM::Matrix(_dimension, bases);
    
    for (unsigned row = 0; row < _dimension; ++row) {
        for (unsigned col = 0; col < bases; ++col) {
            _podBasis(row, col) = (*_U)(row, col);
        }
    }
    
    _indices.push_back(BEM::maxVectorIndex(_podBasis.col(0)));
    for (unsigned col = 1; col < bases; ++col) {
        BEM::Matrix A(col, col);
        BEM::ColVector b(col);
        for (unsigned i = 0; i < col; ++i) {
            b[i] = _podBasis(_indices[i], col);
            for (unsigned j = 0; j < col; ++j) {
                A(i, j) = _podBasis(_indices[i], j);
            }
        }

        BEM::ColVector s(col);
        if (col == 1) {
            s[0] = (b[0]/A(0, 0));
        } else {
            s = A.colPivHouseholderQr().solve(b);
        }
        auto errorVector = _podBasis.block(0, 0, _podBasis.rows(), col)*s - _podBasis.col(col);
        int index = BEM::maxVectorIndex(errorVector);
        if (std::abs(errorVector(index)) < 0.0001) {
            break;
        }
        
        _indices.push_back(index);
    }

    _projMatrix = BEM::Matrix(bases, bases);
    for (unsigned col = 0; col < bases; ++col) {
        for (unsigned row = 0; row < bases; ++row) {
            _projMatrix(row, col) = _podBasis(_indices[row], col);
        }
    }
    msg(5) << "(end)EmpiricalInterpolation::buildInterpolator" << endMsg;
}

/**
 * @brief: build the empirical interpolation with a specified tolerance.
 */
void EmpiricalInterpolation::buildInterpolator(double tolerance)
{
    computeSVD();
    assert(tolerance < 1.);
    msg(1) << "Building interpolator for tolerance: " << tolerance << endMsg;
    std::vector<double> sing(_singVals.size(), 0.0);
    if (not _name.empty()) {
        for (unsigned i = 0; i < sing.size(); ++i) {
            sing[i] = std::log10(std::abs(_singVals[i]/_singVals[0]));
        }
    }
    for (unsigned i = 0; i < sing.size(); ++i) {
        sing[i] = std::pow(std::abs(_singVals[i]), 2);
    }
    double norm = std::accumulate(sing.begin(), sing.end(), 0.0);
    double current = 0.;
    unsigned bases = 0U;
    for (unsigned i = 0U; i < sing.size(); i++) {
        current += sing[i];
        if (1.0 - std::pow(tolerance,2) <= (current/norm)) {
            bases = i + 1;
            break;
        }
    }

    msg(1) << bases << " bases are needed with energy " << std::sqrt(1 - current/norm) << endMsg;
    //assert(bases > 0);
    buildInterpolator(bases);
}

/**
 * @brief: Get the indices necessary for the interpolation.
 */
const std::vector<int> &EmpiricalInterpolation::getInterpolationIndices(void) const
{
    return _indices;
}

/**
 * @brief: Given a vector of values for the necessary indices, get the projection onto the rb.
 */
BEM::ColVector EmpiricalInterpolation::getSolution(const BEM::ColVector &indices) const
{
    assert(indices.size() == _projMatrix.rows());
    if (indices.size() == 0) {
        return BEM::ColVector(0);
    }
    
    return  _projMatrix.colPivHouseholderQr().solve(indices);
}

/**
 * @brief: Overload with std::vector<>
 */
BEM::ColVector EmpiricalInterpolation::getSolution(const std::vector<BEM::Complex> &indices) const
{    
    return  getSolution(BEM::toColVector(indices));
}

/**
 * @brief: Get the resulting vector from the projection onto the rb.
 */
BEM::ColVector EmpiricalInterpolation::interpolate(const BEM::ColVector &indices) const
{
    assert(indices.size() == _projMatrix.rows());
    if (indices.size() == 0) {
        return BEM::ColVector(0);
    }
    
    BEM::ColVector solution = (_projMatrix.colPivHouseholderQr().solve(indices));

    return _podBasis*solution;
}

/**
 * @brief: Overload with std::vector<>
 */
BEM::ColVector EmpiricalInterpolation::interpolate(const std::vector<BEM::Complex> &indices) const
{
    return interpolate(BEM::toColVector(indices));
}
    
