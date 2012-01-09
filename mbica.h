#ifndef MBICA_H
#define MBICA_H

#include<armadillo>

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
        arma::mat count(arma::mat X) {

        }
    };
}
#endif // MBICA_H
