#ifndef POLICIES_H
#define POLICIES_H

#include <armadillo>

namespace mbica {

class NoStabilization {
public:
    NoStabilization(double, double, int) {}
    void operator()(int, const arma::mat &, const arma::mat &) {}
};

class WithStabilization {
public:
    WithStabilization(double &epsilon, double &mu, int &maxIterations)
        : epsilon_(epsilon), mu_(mu), stroke_(0.0), maxIterations_(maxIterations), reducedStep_(false), inStroke_(false)  {}

    void operator()(int iteration, const arma::mat &B, const arma::mat &B_old);

private:
    arma::mat BOlder_;
    double &epsilon_;
    double &mu_;
    double stroke_;
    int maxIterations_;
    bool reducedStep_;
    bool inStroke_;
};
}

#endif // POLICIES_H
