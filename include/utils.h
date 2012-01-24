/**
  * @file
  * @author  Pawel Zubrycki <P.Zubrycki@stud.elka.pw.edu.pl>
  * @author  Stanisław Janikowski <S.A.Janikowski@stud.elka.pw.edu.pl>
  *
  * @section DESCRIPTION
  *
  * File with implementation of algorithms needed to calculate FastICA
  * not provided by Armadillo library.
  **/

#ifndef UTILS_H
#define UTILS_H

#include <armadillo>

namespace mbica {

/**
  * Function to calculate square root of a matrix using diagonalization.
  *
  * @param X Matrix from which will be calculated square root.
  * @return Matrix which will be equal X whan multiplied by itself.
  **/
arma::mat matSqrt(const arma::mat &X);

/**
  * Class implementing primary component analisis (PCA).
  **/
class PCA {
public:
    /**
      * Function that calculaates transformation matrix.
      *
      * @param X Matrix with data on which PCA will be calculated. Samples should be in columns.
      * @param E Output parameter which returns matrix with eigenvectorsin rows, which can be used as transofmration matrix.
      * @param D Output parameter which returns vector with eigenvalues.
      **/
    void operator()(const arma::mat &X, arma::mat &E, arma::vec &D);
};

/**
  * Class implementing whitening transformation
  **/
class Whitening {
public:
    /**
      * Function that calculates whitening/dewhitening matrixes.
      *
      * @param E Matrix with eigenvectors in rows (from PCA).
      * @param D Vector with eigenvalues (from PCA)
      * @param Wh Output parameter which returns whitening matrix.
      * @param dWh Output parameter which returns dewhitening matrix.
      **/
    void operator()(const arma::mat &E, const arma::vec &D,
                    arma::mat &Wh, arma::mat &dWh);
};

/**
  * Function to orthonormalize a matrix.
  *
  * @param A Matrix to orthonormalize.
  * @return Orthonormalized matrix from matrix A.
  **/
arma::mat orth(const arma::mat &A);

/**
  * Function used to remove mean from data.
  *
  * @note Signals have to be in rows, which means samples have to be in columns.
  *
  * @param X Matrix with data.
  * @return Matrix with signals which mean values equal 0.
  **/
arma::mat remmean(const arma::mat &X);
}

#endif // UTILS_H
