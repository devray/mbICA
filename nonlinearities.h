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
            Nonlinearity();

        protected:
            arma::mat g_;
            arma::mat dg_;
        };

        template<int a>
        class Pow: public Nonlinearity{
        public:
            Pow(arma::mat X);
        };

        template<int a, int b>
        class TanH: public Nonlinearity{
            TanH(arma::mat X);
        };

        template<int a, int b>
        class Gauss: public Nonlinearity{
            Gauss(arma::mat X);
        };

        typedef Pow<2> Skew;
    }
}
#endif // NONLINEARITIES_H
