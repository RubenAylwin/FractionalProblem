#ifndef QMC_TEST
#define QMC_TEST

#include <Halton.h>
#include <cmath>
namespace qmcTests{
    static double tolerance = 1E-4;
    BOOST_AUTO_TEST_SUITE(QMC)

    BOOST_AUTO_TEST_SUITE(Halton)

    BOOST_AUTO_TEST_CASE(IntegrationTest1D)
    {
        msg(1) << "start QMC::HALTON::IntgrationTest1D" << endMsg;

        auto fun = [](double t) -> double {
            return std::sqrt(t);
        };
        
        double realResult = 2.*std::sqrt(3.) - 2./3.;
        double error1 = 10;
        double error2 = 10;
        for (int j = 0; j < 10; ++j) {
            double result = 0.0;
            int totPoints = 100*std::pow(2,j);
            for (int i = 0; i < totPoints; ++i) {
                double *point = halton(i, 1);
                result += 2*fun(point[0]*2 + 1)/totPoints;
                delete[] point;
            }
            error2 = std::abs(realResult - result);
            msg(2) << "Error: " << error1/error2 << endMsg;
            BOOST_CHECK(error1/error2 >= 1.7);
            error1 = error2;
        }
        msg(1) << "end QMC::HALTON::IntgrationTest1D" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(IntegrationTest2D)
    {
        msg(1) << "start QMC::HALTON::IntgrationTest2D" << endMsg;

        auto fun = [](double t, double s) -> double {
            return std::sqrt(t)*sqrt(s);
        };
        
        double realResult = std::pow(2.*std::sqrt(3.) - 2./3., 2);
        double error1 = 10;
        double error2 = 10;
        for (int j = 0; j < 10; ++j) {
            double result = 0.0;
            int totPoints = 100*std::pow(2,j);
            for (int i = 0; i < totPoints; ++i) {
                double *point = halton(i, 2);
                result += 2*2*fun(point[0]*2 + 1, point[1]*2 + 1)/totPoints;
                delete[] point;
            }
            error2 = std::abs(realResult - result);
            msg(2) << "Error: " << error1/error2 << endMsg;
            BOOST_CHECK(error1/error2 >= 1.6);
            error1 = error2;
        }
        msg(1) << "end QMC::HALTON::IntgrationTest2D" << endMsg;
    }

    BOOST_AUTO_TEST_CASE(IntegrationTestTrig2D, * boost::unit_test::disabled())
    {
        msg(1) << "start QMC::HALTON::IntgrationTestTrig2D" << endMsg;

        auto fun = [](double t, double s) -> double {
            return std::sin(M_PI*0.5*(t+s));
        };
        
        double realResult = 8./(M_PI*M_PI);
        double result = 0.0;
        int totPoints = 100*std::pow(2,10);
        for (int i = 0; i < totPoints; ++i) {
            double *point = halton(i, 2);
            result += fun(point[0], point[1])/totPoints;
            delete[] point;
        }
        BOOST_CHECK_CLOSE(result, realResult, tolerance);
        msg(1) << "end QMC::HALTON::IntgrationTestTrig2D" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE_END()

}

#endif
