#include <armadillo>
#include <cmath>
using namespace arma;

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "icaseparator.h"

BOOST_AUTO_TEST_SUITE(icaseparator)

BOOST_AUTO_TEST_CASE(separator_multiplying)
{
    mat A = 2 * eye(3, 3);

    mbica::ICASeparator sep(A, A);
    mat B, C, D;
    D << 1 << 2 << 3 << 4 << endr;
    B = repmat(D, 3, 1);
    C << 2 << 4 << 6 << 8 << endr;
    C = repmat(C, 3, 1);

    // Checking if sep(B) == C
    BOOST_CHECK_SMALL(mat(abs(sep(B) - C)).max(), 0.000001);

    // should throw because of bad D size
    BOOST_CHECK_THROW(sep(D), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()

