#ifndef ICASEPARATOR_H
#define ICASEPARATOR_H

#include <armadillo>

namespace mbica {
class ICASeparator
{
public:
    ICASeparator(arma::mat A, arma::mat W)
        : A_(A), W_(W) {A.save("A.mat", arma::arma_ascii); W.save("W.mat", arma::arma_ascii);}

    arma::mat operator()(arma::mat X);

    const arma::mat &getA() const {
        return A_;
    }
    const arma::mat &getW() const {
        return W_;
    }

private:
    arma::mat A_;
    arma::mat W_;
};
}

#endif // ICASEPARATOR_H
