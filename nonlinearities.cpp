#include "nonlinearities.h"

using namespace mbica::nonlinearities;

template<int a>
Pow::Pow(arma::mat X){
    g_ = arma::pow(X, a);
    dg_ = a*arma::pow(X, a-1);
}

template<int a, int b>
TanH::TanH(arma::mat X){
    double c = double(a)/b;
    g_ = arma::tanh(c*X);
    dg_ = c*(1 - pow(g_, 2));
}

template<int a, int b>
Gauss::Gauss(arma::mat X){
    double c = double(a)/b;
    g_ = X * arma::exp(-0.5*c*arma::pow(X, 2));
    dg_ = (1-c*arma::pow(X, 2))*g_;
}



