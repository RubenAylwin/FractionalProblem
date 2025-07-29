#include <Operator.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include <tbb/tbb.h>
#include <string>
#include <Utilities.h>


/**
 * @brief: Constructs operator and initializes mass matrix.
 */
Operator::Operator(const unsigned testSize, const unsigned trialSize) :
    _massMatrix(nullptr),
    _testSize{testSize},
    _trialSize{trialSize}
{
}

/**
 * @brief: Assemble the mass matrix using each childs indexedDuality method.
 * Reads env var to parallelize (or not) the construction.
 */
void Operator::assembleMassMatrix(void)
{
    if (_massMatrix) {
        return;
    }
    _massMatrix.reset(new BEM::Matrix(_testSize, _trialSize));
    if (BEM::getEnv<bool>(std::string("PARALLELIZE_MATRIX_CONSTRUCTION"))) {
        tbb::parallel_for( tbb::blocked_range2d<int>(0, _massMatrix->rows(), 0, _massMatrix->cols()),
                           [this]( const tbb::blocked_range2d<int> &r ) {
                               for(int j=r.rows().begin(), j_end=r.rows().end(); j<j_end; ++j){
                                   for(int k=r.cols().begin(), k_end=r.cols().end(); k<k_end; ++k){
                                       _massMatrix->operator()(j, k) = indexedDuality(j, k);

                                   }

                               }
                           });
    } else {
        for (int i = 0; i < _massMatrix->rows(); ++i) {
            for (int j = 0; j < _massMatrix->cols(); ++j) {
                _massMatrix->operator()(i, j) = indexedDuality(i, j);
            }
        }
    }
}

/**
 * @brief: Returns the mass matrix. Constructs it if it has not been constructed first.
 */
const BEM::Matrix &Operator::getMatrix(void)
{
    if (not _massMatrix) {
        assembleMassMatrix();
    }

    return *_massMatrix;
}
