#include "utils.h"
#include <armadillo>
using namespace arma;

#define BOOST_TEST_MODULE PCA_tests
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(eye_identity)
{
    mat E;
    vec D;
    mat X = eye(10, 10);

    mbica::PCA()(X, E, D);

    BOOST_CHECK(mat(abs(X - eye(10,10))).max() < 0.00001f);
}

BOOST_AUTO_TEST_CASE(negative_eigenvalues)
{
    mat E;
    vec D;
    for(int i=0; i<1000; ++i){
        mat X = randu<mat>(10, 200);
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
