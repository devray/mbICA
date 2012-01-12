#ifndef MBICA_H
#define MBICA_H

#include <armadillo>

//will move to cpp in future
using namespace arma;

namespace mbica {

    class Mbica {
    public:
        Mbica();
    };

    class Pow3 {
    public:
        static arma::mat G(arma::mat X) {
            return arma::pow(X, 3);
        }

        static arma::mat dG(arma::mat X) {
            return 3*arma::pow(X, 2);
        }
    };

    arma::mat matSqrt(const arma:: mat &X) {
        vec D;
        mat E;
        eig_sym(D, E, X);

        return E * sqrt(diagmat(D)) * inv(E);
    }

    template<int mu = 1, class UsedFunc = Pow3>
    class FastICA {
    public:

        void PCA(const arma::mat &X, arma::mat &E, arma::vec &D) {
            mat C = cov(X.t());

            eig_sym(D, E, C);

            int stopIx = D.n_elem, startIx=0;
            for(; stopIx > 0 && D[stopIx-1] == 0; --stopIx);
            for(; startIx < D.n_elem && D[startIx] == 0; ++startIx);

            if(stopIx != D.n_elem || startIx != 0) {
                D = D.subvec(startIx, stopIx-1);
                // TODO: trzeba sprawdzic, czy eigenvectory w wierszach, czy kolumnach
                E = E.cols(startIx, stopIx-1);
            }

            // TODO: sprawdzic, czy nie potrzebna jest transpozycja E - wszak transponujemy na poczatku X
        }

        // E - eigenvec from PCA
        // D - vector with eigenvalues from PCA
        void whitening(const arma::mat &E, const arma::vec &D,
                            arma::mat &Wh, arma::mat &dWh) {
            mat sqrtD = sqrt(diagmat(D));
            Wh = inv(sqrtD) * E.t();
            dWh = E * sqrtD;
        }


        //  tu w jakis posob trzeba umozliwic podanie guess, whitening i dewhitenign matrix
        arma::mat count(arma::mat X, int nIC = -1) {
            double epsilon = 0.1;
            //nIC == -1 means the same size it's now.
            if(nIC == -1) {
                nIC = X.n_rows;
            }

            mat dWh, Wh, E;
            vec D;

            // to możemy zrobić raz i w jakiś sposów podać
            PCA(X, E, D);
            whitening(E ,D, Wh, dWh);

            // Tu whitening (ktory zawiera w sobie PCA)
            mat A = zeros<mat>(X.n_rows, nIC);
            // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
            // B = whiteningMatrix * guess;
            // orth z octave - musimy coś wymyslić
            mat B = arma::orth(randu<mat>(X.n_rows, nIC) - 0.5),
                B_old = arma::zeros<mat>(X.n_rows, nIC);
            int max_iterations = 1000, i;

            // main loop (easiest way - no stabilization, nor fine-tunung)
            for(i=0; i < max_iterations; ++i) {
                B = B * arma::real(matSqrt(arma::inv(B.t() * B)));

                arma::mat minAbsCos = arma::min(arma::abs(arma::diagmat(B.t() * B_old)));
                if( 1 - minAbsCos < epsilon) {
                    break;
                }
                B_old = B;
                double k = 1.0 / X.n_cols;
                arma::mat Y = X.t() * B;
                B = k * (X * UsedFunc::G(Y) - arma::repmat(arma::sum(UsedFunc::dG(Y)), X.n_rows, 1));

            }

            if(i == max_iterations) {
                // return empty A and W
            }
            A = dWh * B;
            W = B.t() * Wh;
            return A, W;
        }
    };
}
#endif // MBICA_H
