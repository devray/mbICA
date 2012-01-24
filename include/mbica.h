#ifndef MBICA_H
#define MBICA_H

#include <armadillo>
#include <boost/parameter.hpp>
#include <iostream>
#include <ctime>

#include "icaseparator.h"
#include "nonlinearities.h"
#include "policies.h"
#include "utils.h"


namespace mbica {

// default values
namespace def {
const double EPSILON = 0.0001;
const int MAX_ITERATIONS = 10000;
const double MU = 1.0;
}

BOOST_PARAMETER_NAME(dWh)
BOOST_PARAMETER_NAME(wh)
BOOST_PARAMETER_NAME(guessMatrix)
BOOST_PARAMETER_NAME(epsilon)
BOOST_PARAMETER_NAME(maxIterations)

// base class for FastICA
class FastICA_impl {
protected:
    template <class ArgumentPack>
    FastICA_impl(ArgumentPack const &args)
        : dWh_(args[_dWh | arma::mat()]),
          Wh_(args[_wh | arma::mat()]),
          guess_(args[_guessMatrix | arma::mat()]),
          epsilon_(args[_epsilon | def::EPSILON]),
          maxIterations_(args[_maxIterations | def::MAX_ITERATIONS])
    {
    }

public:
    void setWhiteningMatrix(arma::mat Wh, arma::mat dWh) {
        Wh_ = Wh;
        dWh_ = dWh;
    }

    void resetWhiteningMatrix() {
        Wh_.reset();
        dWh_.reset();
    }

    const arma::mat &whiteningMatrix() const {
        return Wh_;
    }

    const arma::mat &dewhiteningMatrix() const {
        return dWh_;
    }

    void setEpsilon(double epsilon) {
        epsilon_ = epsilon;
    }

    double epsilon() const {
        return epsilon_;
    }

    void setMaxIterations(int max) {
        maxIterations_ = max;
    }

    int maxIterations() const {
        return maxIterations_;
    }

protected:
    arma::mat dWh_;
    arma::mat Wh_;
    arma::mat guess_;
    double epsilon_;
    int maxIterations_;
};

template<class UsedNonl = nonlinearities::Pow<3>, class Stabilization = NoStabilization >
class FastICA: public FastICA_impl {
public:
    BOOST_PARAMETER_CONSTRUCTOR(
            FastICA, (FastICA_impl), mbica::tag
            , (optional
               (maxIterations, *)
               (epsilon, *)
               (guessMatrix, *)
               (dWh, *)
               (wh, *)
            ))

    ICASeparator operator()(arma::mat X, int nIC = -1, double mu = def::MU) {
        if(Wh_.is_empty() || dWh_.is_empty()) {
            arma::mat E;
            arma::vec D;

            X = remmean(X);
            PCA()(X, E, D);
            std::cout<< D << std::endl;
            Whitening()(E ,D, Wh_, dWh_);
            X = Wh_ * X;
        }

        //nIC == -1 means the same size it's now.
        if(nIC < 0 || unsigned(nIC) > X.n_rows) {
            nIC = X.n_rows;
        }

        // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
        arma::mat B;
        if (!guess_.is_empty())
            B = Wh_ * guess_;
        else
            B = orth(arma::randu<arma::mat>(X.n_rows, nIC) - 0.5);

        arma::mat B_old = arma::zeros<arma::mat>(X.n_rows, nIC);
        Stabilization stabilize(epsilon_, mu, maxIterations_);

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
            if(mu == 1.0) {
                B = k * (X * nl.G()) - k * (arma::repmat(arma::sum(nl.dG()), X.n_rows, 1)) % B;
                //std::cout << B;
            } else {
                arma::mat Beta = sum(Y % nl.G());
                arma::mat D = diagmat(1.0 / (Beta - sum(nl.dG())));
                B += mu * (B * (Y.t() * nl.G() - diagmat(Beta)) * D);
            }
        }
        std::cout << "Number of iters: " << i << std::endl;
        //std::cout << "We got B = " << B << std::endl;

        return ICASeparator(dWh_ * B, B.t() * Wh_);
    }
};
}

#endif // MBICA_H
