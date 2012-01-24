/**
  * @file
  * @author  Pawel Zubrycki <P.Zubrycki@stud.elka.pw.edu.pl>
  * @author  Stanisław Janikowski <S.A.Janikowski@stud.elka.pw.edu.pl>
  *
  * @section DESCRIPTION
  *
  * Class that provides easy way to count independent components.
  **/

#ifndef ICASEPARATOR_H
#define ICASEPARATOR_H

#include <armadillo>

namespace mbica {
/**
  * Class that provides easy way to count independent components.
  **/
class ICASeparator
{
public:
    /**
      * Constructor.
      *
      * @param A Mixing matrix.
      * @param W Unmixing matrix.
      **/
    ICASeparator(arma::mat A, arma::mat W)
        : A_(A), W_(W) {
        A.save("A.mat", arma::arma_ascii);
        W.save("W.mat", arma::arma_ascii);
    }

    /**
      * Method that unmix given signal.
      *
      * @param X Matrix with mixed signal.
      * @return Matrix with unmixed signal.
      **/
    arma::mat operator()(arma::mat X) {
        return W_ * X;
    }

    /**
      * Method that mixes given signals
      *
      * @param S Matrix to be mixed.
      * @return Matrix with mixed data.
      */
    arma::mat mixSignals(arma::mat S) {
        return A_ * S;
    }

    /// Getter of mixing matrix.
    const arma::mat &getA() const {
        return A_;
    }

    /// Getter of unmixing matrix.
    const arma::mat &getW() const {
        return W_;
    }

private:
    arma::mat A_;
    arma::mat W_;
};
}

#endif // ICASEPARATOR_H
