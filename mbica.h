#ifndef MBICA_H
#define MBICA_H

#include <armadillo>

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

    class PCA {
    public:
        static arma::mat process(arma::mat X) {
            // Placeholder to count PCA
            return arma::mat();
        }
    };

    class Whitenig {
        static arma::mat process(arma::mat X) {
            // Placeholder to count PCA and whitening
            return arma::mat();
        }

    };

    template<int mu = 1, class Preprocess = Whitening, class UsedFunc = Pow3>
    class FastICA {
    public:

        arma::mat count(arma::mat X, int nIC = -1) {
            double epsilon = 0.1;
            //nIC == -1 means the same size it's now.
            if(nIC == -1) {
                nIC = X.n_rows;
            }

            X = Preprocess::process(X);
            // Tu whitening (ktory zawiera w sobie PCA)
            arma::mat A = zeros<mat>(X.n_rows, nIC);
            // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
            // B = whiteningMatrix * guess;
            // orth z octave - musimy coś wymyslić
            arma::mat B = arma::orth(randu<mat>(X.n_rows, nIC) - 0.5),
                B_old = arma::zeros<mat>(X.n_rows, nIC);
            int max_iterations = 1000, i;

            // main loop (easiest way)
            for(i=0; i < max_iterations; ++i) {
                // tu trzeba znalezc sposob na znalezienie sqare roota z B
                B = B * arma::real(arma::inv(B.t() * B)^0.5);

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
            A = dewhiteningMatrix * B;
            W = B.t() * whiteningMatrix;
            return A, W;
        }
    };
}
#endif // MBICA_H
