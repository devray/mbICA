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

    cout << sep2.getA() << endl;
    BOOST_CHECK_SMALL(mat(abs(sep2(X) - S)).max(), 0.000001);
}

BOOST_AUTO_TEST_SUITE_END()
