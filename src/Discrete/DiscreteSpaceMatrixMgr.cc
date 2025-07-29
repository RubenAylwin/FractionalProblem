#include <DiscreteSpaceMatrixMgr.h>
#include <DiscreteSpace.h>
#include <Utilities.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include <iostream>


/**
 * @brief: This class is a singleton which will manage all matrices for a multilevel application.
 */
static std::unique_ptr<DiscreteSpaceMatrixMgr_1D> _singleton = nullptr;

/**
 * @brief: Getter
 */
DiscreteSpaceMatrixMgr_1D &DiscreteSpaceMatrixMgr_1D::get(void)
{
    if (not _singleton) {
        _singleton.reset(new DiscreteSpaceMatrixMgr_1D());
    }
    
    return *_singleton;
}

/**
 * @brief: Given a two discrete spaces, it returns the corresponding DOF matrix.
 */
const BEM::Matrix &DiscreteSpaceMatrixMgr_1D::getData(const DiscreteSpace_1D &functionals, const DiscreteSpace_1D &bases)
{
    // Expected to be called from different threads, so lock the data.
    std::lock_guard<std::mutex> lock(_matrixMutex);

    //Get the name of the space. For now we enforce that the spaces must be of the same type, so same name (may come from different discretizations).
    const std::string &functionalsName = functionals.getBaseName();
    // TODO: Change to exception
    assert(functionalsName == bases.getBaseName());
    
    // Construct the corresponding key: Name + size_dof + size_base.
    // Assumption, only one type of discretization is being used OR different discretization strategies have different names.
    std::string key = functionalsName + std::to_string(functionals.getSize()) + "to" +std::to_string(bases.getSize());

    // If the matrix has not been built, build it and save it to the key.
    if (_map.find(key) == _map.end()) {
        _map[key].reset(new BEM::Matrix(functionals.getSize(), bases.getSize()));
        auto & matrix = *(_map[key]);
        if (BEM::getEnv<bool>(std::string("PARALLELIZE_MATRIX_CONSTRUCTION"))) {
            tbb::parallel_for( tbb::blocked_range2d<int>(0, matrix.rows(), 0, matrix.cols()),
                               [&]( const tbb::blocked_range2d<int> &r ) {
                                   for(int j=r.rows().begin(), j_end=r.rows().end(); j<j_end; ++j){
                                       for(int k=r.cols().begin(), k_end=r.cols().end(); k<k_end; ++k){
                                           auto &DoF = functionals.degreeOfFreedom(j);
                                           auto &basis = bases.basisFunction(k);
                                           matrix(j, k) = DoF(basis);
                                       }
                                   }
                               });
        } else {
            for (int i = 0; i < matrix.cols(); ++i) {
                for (int j = 0; j < matrix.rows(); ++j) {
                    auto &DoF = functionals.degreeOfFreedom(j);
                    auto &basis = bases.basisFunction(i);
                    matrix(j, i) = DoF(basis);
                }
            }
        }
    }
    return *_map[key];
}

/**
 * @brief: This method saves to a dictionary all the indices of a given column for which there is non-zero data in the DOF matrix of the given spaces.
 */
const std::vector<int> &DiscreteSpaceMatrixMgr_1D::nonZeroInCol(const DiscreteSpace_1D &from, const DiscreteSpace_1D &to, const unsigned col)
{
    std::lock_guard<std::mutex> lock(_columnMutex);
    const std::string &fromName = from.getBaseName();
    assert(fromName == to.getBaseName());
    std::string key = fromName + std::to_string(from.getSize()) + "to" +std::to_string(to.getSize());
    if (_nonZeroColumns.find(key) == _nonZeroColumns.end()) {
        auto matrix = getData(from, to);
        _nonZeroColumns.insert({key, std::unordered_map<int, std::vector<int>>()});
        for (int col = 0; col < matrix.cols(); ++col) {
            _nonZeroColumns.at(key).insert({col, std::vector<int>()});
            auto column = matrix.col(col);
            for (int i = 0; i < column.size(); ++i) {
                if (column[i] != 0.0) {
                    _nonZeroColumns.at(key).at(col).push_back(i);
                }
            }
        }
    }
    return _nonZeroColumns.at(key).at(col);
}

/**
 * @brief: This method saves to a dictionary all the indices of a given row for which there is non-zero data in the DOF matrix of the given spaces.
 */
const std::vector<int> &DiscreteSpaceMatrixMgr_1D::nonZeroInRow(const DiscreteSpace_1D &from, const DiscreteSpace_1D &to, const unsigned row)
{
    std::lock_guard<std::mutex> lock(_rowMutex);
    const std::string &fromName = from.getBaseName();
    assert(fromName == to.getBaseName());
    std::string key = fromName + std::to_string(from.getSize()) + "to" +std::to_string(to.getSize());
    if (_nonZeroRows.find(key) == _nonZeroRows.end()) {
        auto matrix = getData(from, to);
        _nonZeroRows.insert({key, std::unordered_map<int, std::vector<int>>()});
        for (int row = 0; row < matrix.rows(); ++row) {
            _nonZeroRows.at(key).insert({row, std::vector<int>()});
            auto rowV = matrix.row(row);
            for (int i = 0; i < rowV.size(); ++i) {
                if (rowV[i] != 0.0) {
                    _nonZeroRows.at(key).at(row).push_back(i);
                }
            }
        }
    }
    return _nonZeroRows.at(key).at(row);
}
