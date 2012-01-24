#include <armadillo>
#include <cmath>
using namespace arma;

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "policies.h"

BOOST_AUTO_TEST_SUITE(policies)

BOOST_AUTO_TEST_CASE(stabilization_stroke_check)
{
    double mu = 10.0,
            epsilon = 0.8;
    int max_iterations = 2;
    mat B = 3 * eye(3, 3);
    mat B_old = 0.1 * eye(3, 3);
    mbica::WithStabilization ws(epsilon, mu, max_iterations);

    // setting BOlder to 0.1 * eye(3,3)
    ws(0, B, B_old);

    // checking if stroke will work
    ws(0, B, B_old);
    BOOST_CHECK_SMALL(mu - 5.0, 0.0000001);

    // now mu should go back to 1.0
    ws(0, B, B_old);
    BOOST_CHECK_SMALL(mu - 10.0, 0.0000001);

    B = zeros(3, 3);
}

BOOST_AUTO_TEST_CASE(stabilization_reduceStep_check)
{
    double mu = 10.0,
            epsilon = 0.9;
    int max_iterations = 3;
    mat B = zeros(1, 1);
    mat B_old = zeros(1, 1);
    mbica::WithStabilization ws(epsilon, mu, max_iterations);

    // now mu should be halved
    ws(2, zeros(1, 1), B_old);
    BOOST_CHECK_SMALL(mu - 5.0, 0.0000001);

    // ... and should stay that way even after another late iteration
    ws(2, B, B_old);
    BOOST_CHECK_SMALL(mu - 5.0, 0.0000001);

}

BOOST_AUTO_TEST_SUITE_END()
