#include "icaseparator.h"

using namespace mbica;
using namespace arma;

mat ICASeparator::operator()(mat X) {
    return W_ * X;
}
