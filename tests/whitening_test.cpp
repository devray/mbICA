#include "utils.h"
#include <armadillo>
using namespace arma;

#define BOOST_TEST_MODULE Whitening_tests
#define BOOST_T
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( identity_covariance )
{
    mat E, Wh, dWh;
    vec D;
    for(int i=0; i<100; ++i){
        mat Y,X = randu<mat>(10, 9);
        mbica::PCA()(X,E,D);
        mbica::Whitening()(E,D,Wh,dWh);
        Y = X * Wh;
        BOOST_CHECK(mat(abs(cov(Y)-eye(Y.n_rows,Y.n_cols))).max() < 0.00001f);
    }
}

