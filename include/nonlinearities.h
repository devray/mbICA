/**
  * @file
  * @author  Pawel Zubrycki <P.Zubrycki@stud.elka.pw.edu.pl>
  * @author  Stanisław Janikowski <S.A.Janikowski@stud.elka.pw.edu.pl>
  *
  * @section DESCRIPTION
  *
  * File with preimplemented nonlinearities used in FastICA algorithm.
  *
  * @note New nonlinearites don't have to inherit from Nonlinearity class.
  *       They just have to provide G() and dG() methods.
  **/

#ifndef NONLINEARITIES_H
#define NONLINEARITIES_H

#include<armadillo>

namespace mbica {
/// namespace for available nonlinearities
namespace nonlinearities {

/**
  * Base class for preimplemented nonlinearities.
  */
class Nonlinearity {
public:
    /// Method returning value of function.
    arma::mat G() { return g_; }
    /// Method returning value of derivative of function.
    arma::mat dG() { return dg_; }

protected:
    /// Protected constructor to prevent from creating Nonlinearity objects.
    Nonlinearity(){}

protected:
    arma::mat g_;
    arma::mat dg_;
};

/**
  * Class implementing power function on matrix elements.
  *
  * @param a Exponent.
  **/
template<int a=3>
class Pow: public Nonlinearity{
public:
    Pow(arma::mat X){
        g_ = arma::pow(X, a);
        dg_ = a*arma::pow(X, a-1);
    }
};


/**
  * Class implementing tanh nonlinearity.
  *
  * @note Class uses fraction a/b as a parameter.
  *
  * @param a Numerator of parameter.
  * @param b Denominator of parameter.
  **/
template<int a=1, int b=1>
class TanH: public Nonlinearity{
public:
    TanH(arma::mat X){
        double c = double(a)/b;
        g_ = arma::tanh(c * X);
        dg_ = c * (1 - pow(g_, 2));
    }
};

/**
  * Class implementing gauss nonlinearity.
  *
  * @note Class uses fraction a/b as a parameter.
  *
  * @param a Numerator of parameter.
  * @param b Denominator of parameter.
  **/
template<int a=1, int b=1>
class Gauss: public Nonlinearity{
public:
    Gauss(arma::mat X){
        double c = double(a)/b;
        arma::mat ex = arma::exp(-0.5 * c * arma::pow(X, 2));
        g_ = X % ex;
        dg_ = (1 - c * arma::pow(X, 2)) % ex;
    }
};

/// Typedef to provide compatibility with Octave packet.
typedef Pow<2> Skew;
}
}
#endif // NONLINEARITIES_H
