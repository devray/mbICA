#include <armadillo>
#include <cmath>
using namespace arma;

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "nonlinearities.h"

BOOST_AUTO_TEST_SUITE(nonlinearities)

BOOST_AUTO_TEST_CASE(pow3)
{
    mat A, G, dG;
    A  << 1 << 2  << 3  << endr;
    G  << 1 << 8  << 27 << endr;
    dG << 3 << 12 << 27 << endr;

    mbica::nonlinearities::Pow<3> nl(A);

    BOOST_CHECK_SMALL(mat(abs(nl.G() - G)).max(), 0.000001);
    BOOST_CHECK_SMALL(mat(abs(nl.dG() - dG)).max(), 0.000001);
}

BOOST_AUTO_TEST_CASE(tanH)
{
    mat A, G, dG;
    A  << 1 << 2 << 3  << endr;
    G  << 0.321513 << 0.582783 << 0.761594 << endr;
    dG << 0.298877 << 0.220121 << 0.139991 << endr;

    mbica::nonlinearities::TanH<1,3> nl(A);

    BOOST_CHECK_SMALL(mat(abs(nl.G() - G)).max(), 0.000001);
    BOOST_CHECK_SMALL(mat(abs(nl.dG() - dG)).max(), 0.000001);
}

BOOST_AUTO_TEST_CASE(gauss)
{
    mat A, G, dG;
    A  << 1 << 2  << 3  << endr;
    G  << 0.846482 << 1.02683 << 0.66939 << endr;
    dG << 0.564321 << -0.171139 << -0.44626 << endr;

    mbica::nonlinearities::Gauss<1,3> nl(A);

    BOOST_CHECK_SMALL(mat(abs(nl.G() - G)).max(), 0.00001);
    BOOST_CHECK_SMALL(mat(abs(nl.dG() - dG)).max(), 0.00001);
}


BOOST_AUTO_TEST_SUITE_END()
