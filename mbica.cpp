#include "mbica.h"

using namespace arma;

namespace mbica {

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
    //whiteningMatrixIsSet_ = (dWh() != 0 && Wh != 0);
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
    stroke_ = 0.0;
    reducedStep_ = false;
    guessProvided_ = false;
    srand(time(NULL));
}

void FastICA_impl::stabilize(int iteration, arma::mat B, arma::mat B_older ){
    if (stroke_!=0.0)
        mu_=stroke_;
    else{
        double minAbsCos2 = arma::min(arma::min(arma::abs(arma::diagvec(B.t() * B_older))));
        if (1 - minAbsCos2 < epsilon_){
            stroke_ = mu_;
            mu_ /= 2;
        }
        else if (!reducedStep_ && (iteration > maxIterations_/2)){
            reducedStep_=true;
            mu_ /= 2;
        }
    }

}
}
