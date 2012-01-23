#include "policies.h"

using namespace arma;
using namespace mbica;

void WithStabilization::operator ()(int iteration, const arma::mat &B, const arma::mat &B_old) {
    if (stroke_ != 0.0)
        mu_ = stroke_;
    else {
        if (1.0 - mat(abs(diagvec(B.t() * BOlder_))).min() < epsilon_) {
            stroke_ = mu_;
            mu_ /= 2;
        } else if (!reducedStep_ && (iteration > maxIterations_ / 2)) {
            reducedStep_ = true;
            mu_ /= 2;
        }
    }
    BOlder_ = B_old;
}

