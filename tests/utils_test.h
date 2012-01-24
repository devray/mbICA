#include "utils.h"

// matSqrt test
BOOST_AUTO_TEST_CASE(sqrt_square_equal)
{
    mat X;
    X << 1 << 2 << endr <<
         3 << 4 << endr <<
         6 << 6 << endr;

    mat C = cov(X);
    mat Csqrt = mbica::matSqrt(C);

    BOOST_CHECK_SMALL(mat(abs(Csqrt * Csqrt - C)).max(), 0.000001);
}

// PCA tests
BOOST_AUTO_TEST_CASE(eye_identity)
{
    mat E;
    vec D;
    mat X(3, 6);
    X.cols(0,2) = eye(3, 3);
    X.cols(3,5) = -1 * eye(3, 3);

    mbica::PCA()(X, E, D);

    mat Y = E * X;

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

// Whitening tests
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

