#include <armadillo>
#include <cmath>
using namespace arma;

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "mbica.h"
#include "icaseparator.h"

BOOST_AUTO_TEST_SUITE(fastica)

BOOST_AUTO_TEST_CASE(unmixing)
{
    mat A, W;
    A << 1 << 1 << 0 << endr <<
         0 << 1 << 1 << endr <<
         1 << 0 << 0 << endr;

    W = inv(A);

    mbica::ICASeparator sep(A, W);

    mat S = randu<mat>(3, 1000) - 0.5;
    mat X = sep.mixSignals(S);

    mbica::ICASeparator sep2 =
            mbica::FastICA<mbica::nonlinearities::Gauss<1, 1>, mbica::WithStabilization>()(X);

    cout << "Real A = \n" << (A >= 0.5)  << endl;
    // Ugly way of getting found mixing matrix to former mixing matrix.
    cout << "Found A = \n" << (abs(sep2.getA()) >= 0.2) << endl;
    // K is now a shufled eye matrix
    mat K = inv(A) * (abs(sep2.getA()) >= 0.2);

    // K should be shuffled eye if
    BOOST_CHECK_SMALL(mat(abs(sep2.getA() * sep2.getW() - eye(3,3))).max(), 0.01);

    //Checking if matrixes are similar
    BOOST_WARN_SMALL(sum(sum(K % K)) - 3.0 , 0.00001);
    BOOST_WARN_SMALL(sum(sum(K)) - 3.0 , 0.00001);
}

BOOST_AUTO_TEST_SUITE_END()
