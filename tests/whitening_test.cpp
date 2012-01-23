#include "utils.h"
#include <armadillo>
using namespace arma;

#define BOOST_TEST_MODULE Whitening_tests
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( dontknowyet )
{
    mat E;
    vec D;
    for(int i=0; i<100; ++i){
        mat X = randu<mat>(10, 9);
//        mbica::PCA()(X,E,D);

//        BOOST_CHECK(D.min() >= 0);
    }
}

