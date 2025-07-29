#ifndef UTILS_TEST
#define UTILS_TEST

#include <Utilities.h>
#include <MyTypes.h>

namespace UtilsTest{
    static double tolerance = 1E-3;
    BOOST_AUTO_TEST_SUITE(Utils)

    
    BOOST_AUTO_TEST_CASE(HankelTest) {
        msg(1) << "start Utils::HankelTest" << endMsg;
        double z = 100;
        for (int i=1; i < 100; ++i){
            BOOST_CHECK_CLOSE(BEM::hankel0_1(z).real(), BEM::hankel0_1_orig(z).real(), tolerance);
            BOOST_CHECK_CLOSE(BEM::hankel0_1(z).imag(), BEM::hankel0_1_orig(z).imag(), tolerance);
            --z;
        }
        msg(1) << "end Utils::HankelTest" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(MatrixProductReal) {
        msg(1) << "start Utils::MatrixProductReal" << endMsg;
        BEM::Matrix matrix1(2, 3);
        BEM::Matrix matrix2(3, 2);
        matrix1 << 0.0, 1.0, 0.3,
            1.0, -2.0, 3.0;
        matrix2 << 1., 1.,
            0., -0.5,
            3., 0.;
        BEM::Matrix result = BEM::product(matrix1, matrix2);
        BOOST_CHECK_CLOSE(result(0,0).real(), 0.9, 1E-10);
        BOOST_CHECK_CLOSE(result(0,1).real(), -0.5, 1E-10);
        BOOST_CHECK_CLOSE(result(1,0).real(), 10., 1E-10);
        BOOST_CHECK_CLOSE(result(1,1).real(), 2., 1E-10);
        msg(1) << "end Utils::MatrixProductReal" << endMsg;
    }


    BOOST_AUTO_TEST_CASE(MatrixProductRealComplex) {
        msg(1) << "start Utils::MatrixProductRealComplex" << endMsg;
        BEM::Matrix matrix1(2, 3);
        BEM::Matrix matrix2(3, 2);
        matrix1 << 0.0 + BEM::I, 1.0 - 0.3*BEM::I, 0.3,
            1.0, -2.0 + 0.5*BEM::I, 3.0 - 2.*BEM::I;
        matrix2 << 1., 1. + 0.5*BEM::I,
            BEM::I, -0.5,
            3. - 2.*BEM::I, 0.;
        BEM::Matrix result = BEM::product(matrix1, matrix2);
        BEM::Matrix eigenResult = matrix1*matrix2;
        BOOST_CHECK_CLOSE(result(0,0).real(), eigenResult(0,0).real(), 1E-10);
        BOOST_CHECK_CLOSE(result(0,1).real(), eigenResult(0,1).real(), 1E-10);
        BOOST_CHECK_CLOSE(result(1,0).real(), eigenResult(1,0).real(), 1E-10);
        BOOST_CHECK_CLOSE(result(1,1).real(), eigenResult(1,1).real(), 1E-10);

        BOOST_CHECK_CLOSE(result(0,0).imag(), eigenResult(0,0).imag(), 1E-10);
        BOOST_CHECK_CLOSE(result(0,1).imag(), eigenResult(0,1).imag(), 1E-10);
        BOOST_CHECK_CLOSE(result(1,0).imag(), eigenResult(1,0).imag(), 1E-10);
        BOOST_CHECK_CLOSE(result(1,1).imag(), eigenResult(1,1).imag(), 1E-10);
        msg(1) << "end Utils::MatrixProductRealComplex" << endMsg;
    }

    
    BOOST_AUTO_TEST_SUITE_END()
}

#endif
