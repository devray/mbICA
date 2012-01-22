#ifndef MBICA_H
#define MBICA_H

#include <armadillo>
#include <iostream>
#include "nonlinearities.h"
#include "icaseparator.h"
#include "utils.h"

//will move to cpp in future
using namespace arma;

namespace mbica {

const double DEF_EPSILON = 0.0001;
const int DEF_MAX_ITER = 10000;
const double DEF_MU = 1.0;

class FastICA_impl {
protected:
    FastICA_impl(double epsilon = DEF_EPSILON,
                 int maxIterations = DEF_MAX_ITER,
                 double mu = DEF_MU);
    FastICA_impl(arma::mat Wh, mat dWh);
    FastICA_impl(arma::mat guess);

public:
    void setWhiteningMatrix(arma::mat Wh, mat dWh);

    void unsetWhiteningMatrix() {
        whiteningMatrixIsSet_ = false;
    }

    void setEpsilon(double epsilon) {
        epsilon_ = epsilon;
    }

    void setMaxIterations(int max) {
        maxIterations_ = max;
    }

protected:
    void init(double epsilon = DEF_EPSILON,
              int maxIterations = DEF_MAX_ITER,
              double mu = DEF_MU);

    void stabilize(int iteration, arma::mat B, arma::mat B_older );

protected:
    bool whiteningMatrixIsSet_;
    bool stabilizationEnabled_;
    bool guessProvided_;
    mat dWh_;
    mat Wh_;
    mat guess_;
    double epsilon_;
    int maxIterations_;
    double mu_;
    double stroke_;
    bool reducedStep_;
};

template<class UsedNonl = nonlinearities::Pow<3> >
class FastICA: public FastICA_impl {
public:
    FastICA(double epsilon = DEF_EPSILON, int maxIterations = DEF_MAX_ITER, double mu = DEF_MU)
        : FastICA_impl(epsilon, maxIterations, mu) {}

    FastICA(arma::mat Wh, arma::mat dWh)
        : FastICA_impl(Wh, dWh) {}

    FastICA(arma::mat guess)
        : FastICA_impl(guess) {}

    ICASeparator operator()(arma::mat X, int nIC = -1) {
        //nIC == -1 means the same size it's now.
        if(nIC == -1) {
            nIC = X.n_rows;
        }

        X = remmean(X);

        if(!whiteningMatrixIsSet_) {
            mat E;
            vec D;

            PCA()(X, E, D);
            Whitening()(E ,D, Wh_, dWh_);
        }

        // Tu whitening (ktory zawiera w sobie PCA)
        mat A = zeros<arma::mat>(X.n_rows, nIC);
        // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
        mat B;
        if (guessProvided_)
            B = Wh_ * guess_;
        else
            B = orth(randu<arma::mat>(X.n_rows, nIC) - 0.5);

        mat B_old = arma::zeros<arma::mat>(X.n_rows, nIC);
        mat B_older = arma::zeros<arma::mat>(X.n_rows, nIC);
        int i;

        // main loop (with stabilization, but no fine-tunung)
        for(i = 0; i < maxIterations_; ++i) {
            B = B * matSqrt(arma::inv(B.t() * B));

            double minAbsCos = arma::min(arma::min(arma::abs(arma::diagvec(B.t() * B_old))));
            std::cout << "Step: " << i << ", estiarma::mate: " << 1-minAbsCos << std::endl;
            if( 1 - minAbsCos < epsilon_) {
                break;
            }
            else if(stabilizationEnabled_)
                stabilize(i, B, B_older);

            B_older = B_old;
            B_old = B;
            double k = 1.0 / X.n_cols;
            arma::mat Y = X.t() * B;
            UsedNonl nl(Y);
            if(mu_==1.0)
                B = k * X * nl.G() - k * arma::repmat(arma::sum(nl.dG()), X.n_rows, 1) % B;
            else{
                arma::mat Beta = sum(Y % nl.G());
                arma::mat D = diagmat(1/(Beta-sum(nl.dG())));
                B = B + mu_ * B * (Y.t() * nl.G() - diagmat(Beta))* D;
            }


        }
        std::cout << "Number of iters: " << i << std::endl;
        if(i == maxIterations_) {
            // return empty A and W
        }

        return ICASeparator(dWh_ * B, B.t() * Wh_);
    }
};
}

#endif // MBICA_H
