#ifndef NONLINEARITIES_H
#define NONLINEARITIES_H

#include<armadillo>

namespace mbica {
namespace nonlinearities {

class Nonlinearity {
public:
    arma::mat G() { return g_; }
    arma::mat dG() { return dg_; }

protected:
    Nonlinearity(){}

protected:
    arma::mat g_;
    arma::mat dg_;
};

template<int a>
class Pow: public Nonlinearity{
public:
    Pow(arma::mat X){
        g_ = arma::pow(X, a);
        dg_ = a*arma::pow(X, a-1);
    }
};

template<int a, int b>
class TanH: public Nonlinearity{
    TanH(arma::mat X){
        double c = double(a)/b;
        g_ = arma::tanh(c*X);
        dg_ = c*(1 - pow(g_, 2));
    }
};

template<int a, int b>
class Gauss: public Nonlinearity{
    Gauss(arma::mat X){
        double c = double(a)/b;
        g_ = X * arma::exp(-0.5*c*arma::pow(X, 2));
        dg_ = (1-c*arma::pow(X, 2))*g_;
    }

};

typedef Pow<2> Skew;
}
}
#endif // NONLINEARITIES_H
