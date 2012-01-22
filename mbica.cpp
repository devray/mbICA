#include "mbica.h"

using namespace arma;
using namespace mbica;

FastICA_impl::FastICA_impl(double epsilon, int maxIterations, double mu)
    : whiteningMatrixIsSet_(false) {
    init(epsilon, maxIterations, mu);
}

FastICA_impl::FastICA_impl(arma::mat Wh, mat dWh) {
    init();
    setWhiteningMatrix(Wh, dWh);
}

FastICA_impl::FastICA_impl(arma::mat guess) {
    init();
    guess_=guess;
    guessProvided_=true;
}

void FastICA_impl::setWhiteningMatrix(arma::mat Wh, mat dWh) {
    whiteningMatrixIsSet_ = true;
    if(whiteningMatrixIsSet_){
        dWh_ = dWh;
        Wh_ = Wh;
    }
}

void FastICA_impl::init(double epsilon, int maxIterations, double mu) {
    epsilon_ = epsilon;
    maxIterations_ = maxIterations;
    mu_ = mu;
    guessProvided_ = false;
    srand(time(NULL));
}

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
