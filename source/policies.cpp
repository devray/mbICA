#include "policies.h"

using namespace arma;
using namespace mbica;

void WithStabilization::operator ()(int iteration, const arma::mat &B, const arma::mat &B_old) {
    if ( BOlder_.is_empty() ) BOlder_.copy_size(B_old);
    if (inStroke_){
        mu_ = stroke_;
        inStroke_=false;
    } else {
        if (1.0 - mat(abs(diagvec(B.t() * BOlder_))).min() < epsilon_) {
            stroke_ = mu_;
            mu_ /= 2;
            inStroke_=true;
            std::cout << "Stroke!" << endl;
        } else if (!reducedStep_ && (iteration > maxIterations_ / 2)) {
            reducedStep_ = true;
            mu_ /= 2;
            std::cout << "Reducing!" << endl;
        }
    }
    BOlder_ = B_old;
}

