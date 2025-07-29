#ifndef EMPIRICAL_TEST
#define EMPIRICAL_TEST
#include <MyTypes.h>
#include <DiscreteSpace.h>
#include <DiscreteSpaceMatrixMgr.h>
#include <EmpiricalInterpolation.h>
#include <vector>
#include <chrono>
#include <PECGratingReduced.h>
#include <unordered_set>
#include <progressbar.hpp>
#include <tbb/parallel_for.h>


namespace EmpiricalTest {
    static double tolerance = 1E-10;
    BOOST_AUTO_TEST_SUITE(Empirical)

    BOOST_AUTO_TEST_SUITE(General)

    BOOST_AUTO_TEST_CASE(BasicInterpolation) {
        msg(1) << "start Empirical::General::BasicInterpolation" << endMsg;
        EmpiricalInterpolation emp(2);
        emp.loadTrainingData(std::vector<BEM::Complex>{1.0, 0.0});
        emp.loadTrainingData(std::vector<BEM::Complex>{1.0, 30.0});
        emp.buildInterpolator(0.0000000001);
        std::vector<BEM::Complex> data{2., 3.};
        auto indices = emp.getInterpolationIndices();
        std::vector<BEM::Complex> interp{data[indices[0]], data[indices[1]]};
        auto solution = emp.interpolate(std::move(interp));
        BOOST_CHECK_CLOSE(std::real(solution[0]), 2., tolerance);
        BOOST_CHECK_CLOSE(std::real(solution[1]), 3., tolerance);
        msg(1) << "end Empirical::General::BasicInterpolation" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(BasicInterpolationComplex) {
        msg(1) << "start Empirical::General::BasicInterpolationComplex" << endMsg;
        EmpiricalInterpolation emp(2);
        emp.loadTrainingData(std::vector<BEM::Complex>{BEM::Complex(0, 1.0), 0.0});
        emp.loadTrainingData(std::vector<BEM::Complex>{1.0, 30.0});
        emp.buildInterpolator(0.0000000001);
        std::vector<BEM::Complex> data{2., 3.};
        auto indices = emp.getInterpolationIndices();
        std::vector<BEM::Complex> interp{data[indices[0]], data[indices[1]]};
        auto solution = emp.interpolate(std::move(interp));
        BOOST_CHECK_CLOSE(std::real(solution[0]), 2., tolerance);
        BOOST_CHECK_CLOSE(std::real(solution[1]), 3., tolerance);
        msg(1) << "end Empirical::General::BasicInterpolationComplex" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(BasicInterpolationMoreDim) {
        msg(1) << "start Empirical::General::BasicInterpolationMoreDim" << endMsg;
        EmpiricalInterpolation emp(3);
        emp.loadTrainingData(std::vector<BEM::Complex>{1.0, 0.0, 2.0});
        emp.loadTrainingData(std::vector<BEM::Complex>{1.0, 30.0, -1.0});
        emp.buildInterpolator(0.0000000001);
        std::vector<BEM::Complex> data{2., 3., 0.};
        auto indices = emp.getInterpolationIndices();
        std::vector<BEM::Complex> interp{data[indices[0]], data[indices[1]]};
        auto solution = emp.interpolate(std::move(interp));
        BOOST_CHECK_CLOSE(std::real(solution[1]), 3., tolerance);
        BOOST_CHECK(std::abs(std::real(solution[2])) < 1e-10);
        msg(1) << "end Empirical::General::BasicInterpolationMoreDim" << endMsg;
    }
    BOOST_AUTO_TEST_SUITE_END()    
    
    BOOST_AUTO_TEST_SUITE_END()
}


#endif
