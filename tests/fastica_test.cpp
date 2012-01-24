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

    mat Ap = sep2.getA();
    // normalization
    A = A / repmat(sqrt(sum(pow(A, 2))), A.n_rows, 1);
    Ap = Ap / repmat(sqrt(sum(pow(Ap, 2))), Ap.n_rows, 1);
    cout << "Real A = \n" << A  << endl;
    // Ugly way of getting found mixing matrix to former mixing matrix.
    cout << "Found A = \n" << Ap << endl; //(abs(sep2.getA()) >= 0.2) << endl;
    // K is now a shufled eye matrix
    mat K = (abs(inv(A) * Ap));

    // K should be shuffled eye if just columns are shuffled in both matrixes
    BOOST_CHECK_SMALL(mat(abs(sep2.getA() * sep2.getW() - eye(3,3))).max(), 0.01);

    //Checking if matrixes are similar
    BOOST_CHECK_SMALL(sum(sum(pow(K, 2))) - 3.0 , 0.1);
}

BOOST_AUTO_TEST_SUITE_END()
