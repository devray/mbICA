#include "utils.h"
#include <armadillo>
#include <cmath>
using namespace arma;

#define BOOST_TEST_MODULE PCA_tests
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(eye_identity)
{
    mat E;
    vec D;
    mat X(3, 6);
    X.cols(0,2) = eye(3, 3);
    X.cols(3,5) = -1 * eye(3, 3);

    mbica::PCA()(X, E, D);

    mat Y = E.t() * X;

    BOOST_CHECK_SMALL(mat(abs(Y - X)).max(), 0.000001);
}

BOOST_AUTO_TEST_CASE(negative_eigenvalues)
{
    mat E;
    vec D;
    for(int i=0; i<100; ++i){
        mat X = randu<mat>(10, 9);
        mbica::PCA()(X,E,D);
        BOOST_CHECK(D.min() >= 0);
    }
}

BOOST_AUTO_TEST_CASE(dimension_reduction)
{
    mat E;
    vec D, row;
    row << 1 << 2 << 3 << 4 << 5 << endr;
    mat X = repmat(row, 3, 1);
    mbica::PCA()(X, E, D);

    BOOST_CHECK(D.n_elem == 1);
}
