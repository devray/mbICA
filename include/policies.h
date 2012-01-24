/**
  * @file
  * @author  Pawel Zubrycki <P.Zubrycki@stud.elka.pw.edu.pl>
  * @author  Stanisław Janikowski <S.A.Janikowski@stud.elka.pw.edu.pl>
  *
  * @section DESCRIPTION
  *
  * File with policies that can be applied to FastICA template.
  **/

#ifndef POLICIES_H
#define POLICIES_H

#include <armadillo>

namespace mbica {

/**
  * Empty class that doesn't provide stabilization.
  **/
class NoStabilization {
public:
    NoStabilization(double, double, int) {}
    void operator()(int, const arma::mat &, const arma::mat &) {}
};

/**
  * Class that provides stabilization.
  **/
class WithStabilization {
public:
    /**
      * Constructor
      *
      * @param epsilon Reference to epsilon stored in FastICA class.
      * @param mu Reference to current mu value from FastICA::operator() method.
      * @param maxIterations Reference to value stored in FastICA class.
      **/
    WithStabilization(double &epsilon, double &mu, int &maxIterations)
        : epsilon_(epsilon), mu_(mu), stroke_(0.0), maxIterations_(maxIterations), reducedStep_(false), inStroke_(false)  {}

    /**
      * Method that stabilize calcualtions.
      *
      * @param iteration Current iteration number in main loop of algorithm.
      * @param B Current estimation.
      * @param B_old Previous estimation.
      **/
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
