#ifndef UTILS_H
#define UTILS_H

#include <armadillo>

namespace mbica {
    arma::mat matSqrt(const arma:: mat &X);

    class PCA {
    public:
        void operator()(const arma::mat &X, arma::mat &E, arma::vec &D);
    };

    class Whitening {
    public:
        // E - eigenvec from PCA
        // D - vector with eigenvalues from PCA
        void operator()(const arma::mat &E, const arma::vec &D,
                        arma::mat &Wh, arma::mat &dWh);
    };

    arma::mat orth(const arma::mat &A);

    arma::mat remmean(const arma::mat &X);
}

#endif // UTILS_H
