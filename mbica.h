#ifndef MBICA_H
#define MBICA_H

#include <armadillo>
#include <iostream>
#include "nonlinearities.h"
#include "icaseparator.h"
#include "utils.h"
#include <boost/parameter.hpp>

namespace mbica {

const double DEF_EPSILON = 0.0001;
const int DEF_MAX_ITER = 10000;
const double DEF_MU = 1.0;

BOOST_PARAMETER_NAME(dWh)
BOOST_PARAMETER_NAME(wh)
BOOST_PARAMETER_NAME(guessMatrix)
BOOST_PARAMETER_NAME(epsilon)
BOOST_PARAMETER_NAME(maxIterations)
BOOST_PARAMETER_NAME(mu)

class FastICA_impl {
protected:
    template <class ArgumentPack>
    FastICA_impl(ArgumentPack const &args)
        : dWh_(args[_dWh | arma::mat()]),
          Wh_(args[_wh | arma::mat()]),
          guess_(args[_guessMatrix | arma::mat()]),
          epsilon_(args[_epsilon | DEF_EPSILON]),
          maxIterations_(args[_maxIterations | DEF_MAX_ITER]),
          mu_(args[_mu | DEF_MU])
    {
    }

public:
    void setWhiteningMatrix(arma::mat Wh, arma::mat dWh) {
        Wh_ = Wh;
        dWh_ = dWh;
    }

    void unsetWhiteningMatrix() {
        Wh_.reset();
        dWh_.reset();
    }

    void setEpsilon(double epsilon) {
        epsilon_ = epsilon;
    }

    void setMaxIterations(int max) {
        maxIterations_ = max;
    }

protected:
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
    BOOST_PARAMETER_CONSTRUCTOR(
            FastICA, (FastICA_impl), mbica::tag
            , (optional
               (mu, *)
               (maxIterations, *)
               (epsilon, *)
               (guessMatrix, *)
               (dWh, *)
               (wh, *)
            ))

    ICASeparator operator()(arma::mat X, int nIC = -1) {
        //nIC == -1 means the same size it's now.
        if(nIC == -1) {
            nIC = X.n_rows;
        }

        X = remmean(X);

        if(Wh_.is_empty() || dWh_.is_empty()) {
            arma::mat E;
            arma::vec D;

            PCA()(X, E, D);
            Whitening()(E ,D, Wh_, dWh_);

        }
        X = Wh_ * X;

        // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
        arma::mat B;
        if (!guess_.is_empty())
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
