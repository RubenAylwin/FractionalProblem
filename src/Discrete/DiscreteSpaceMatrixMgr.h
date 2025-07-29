#ifndef DISCRETE_SPACE_MATRIX_MANAGER
#define DISCRETE_SPACE_MATRIX_MANAGER

#include <MyTypes.h>
#include <memory>
#include <unordered_map>
#include <string>
#include <mutex>

////////////////////////////////////////////////////////////////////////
// Class for Matrix Management for MultiLevel Empirical Interpolation //
////////////////////////////////////////////////////////////////////////

class DiscreteSpace_1D;

class DiscreteSpaceMatrixMgr_1D {
public:
    static DiscreteSpaceMatrixMgr_1D &get(void);
    const BEM::Matrix &getData(const DiscreteSpace_1D &from, const DiscreteSpace_1D &to);
    const std::vector<int> &nonZeroInCol(const DiscreteSpace_1D &from, const DiscreteSpace_1D &to, const unsigned col);
    const std::vector<int> &nonZeroInRow(const DiscreteSpace_1D &from, const DiscreteSpace_1D &to, const unsigned row);
private:
    DiscreteSpaceMatrixMgr_1D(void) = default;
    std::unordered_map<std::string, std::unique_ptr<BEM::Matrix>> _map;
    std::mutex _matrixMutex;
    std::mutex _columnMutex;
    std::mutex _rowMutex;
    std::unordered_map<std::string, std::unordered_map<int, std::vector<int>>> _nonZeroColumns;
    std::unordered_map<std::string, std::unordered_map<int, std::vector<int>>> _nonZeroRows;
};

#endif
