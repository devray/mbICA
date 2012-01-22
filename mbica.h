#ifndef MBICA_H
#define MBICA_H

#include <armadillo>
#include <iostream>
#include "nonlinearities.h"
#include "icaseparator.h"
#include "utils.h"

namespace mbica {

const double DEF_EPSILON = 0.0001;
const int DEF_MAX_ITER = 10000;
const double DEF_MU = 1.0;

class FastICA_impl {
protected:
    FastICA_impl(double epsilon = DEF_EPSILON,
                 int maxIterations = DEF_MAX_ITER,
                 double mu = DEF_MU);
    FastICA_impl(arma::mat Wh, arma::mat dWh);
    FastICA_impl(arma::mat guess);

public:
    void setWhiteningMatrix(arma::mat Wh, arma::mat dWh);

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

    void stabilize(int iteration, const arma::mat &B, const arma::mat &B_older );

protected:
    bool whiteningMatrixIsSet_;
    bool guessProvided_;
    arma::mat dWh_;
    arma::mat Wh_;
    arma::mat guess_;
    double epsilon_;
    int maxIterations_;
    double mu_;
};

class NoStabilization {
public:
    NoStabilization(double, double, int) {}
    void operator()(int, const arma::mat &, const arma::mat &) {}
};

class WithStabilization {
public:
    WithStabilization(double &epsilon, double &mu, int &maxIterations)
        : epsilon_(epsilon), mu_(mu), stroke_(0.0), maxIterations_(maxIterations), reducedStep_(false)  {}

    void operator()(int iteration, const arma::mat &B, const arma::mat &B_old);

private:
    arma::mat BOlder_;
    double &epsilon_;
    double &mu_;
    double stroke_;
    int maxIterations_;
    bool reducedStep_;
};

template<class UsedNonl = nonlinearities::Pow<3>, class Stabilization = NoStabilization >
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
            arma::mat E;
            arma::vec D;

            PCA()(X, E, D);
            Whitening()(E ,D, Wh_, dWh_);
        }

        // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
        arma::mat B;
        if (guessProvided_)
            B = Wh_ * guess_;
        else
            B = orth(arma::randu<arma::mat>(X.n_rows, nIC) - 0.5);

        arma::mat B_old = arma::zeros<arma::mat>(X.n_rows, nIC);
        Stabilization stabilize(epsilon_, mu_, maxIterations_);

        double k = 1.0 / X.n_cols;
        int i;
        // main loop (with stabilization, but no fine-tunung)
        for(i = 0; i < maxIterations_; ++i) {
            B = B * matSqrt(arma::inv(B.t() * B));

            double minAbsCos = 1.0 - arma::mat(arma::abs(arma::diagvec(B.t() * B_old))).min();
            std::cout << "Step: " << i << ", estimate: " << minAbsCos << std::endl;
            if(minAbsCos < epsilon_)
                break;

            stabilize(i, B, B_old);

            B_old = B;

            arma::mat Y = X.t() * B;
            UsedNonl nl(Y);
            if(mu_ == 1.0) {
                B = k * (X * nl.G()) - k * (arma::repmat(arma::sum(nl.dG()), X.n_rows, 1) % B);
            } else {
                arma::mat Beta = sum(Y % nl.G());
                arma::mat D = diagmat(1.0 / (Beta - sum(nl.dG())));
                B += mu_ * (B * (Y.t() * nl.G() - diagmat(Beta)) * D);
            }
        }
        std::cout << "Number of iters: " << i << std::endl;
        std::cout << "We got B = " << B << std::endl;

        return ICASeparator(dWh_ * B, B.t() * Wh_);
    }
};
}

#endif // MBICA_H
