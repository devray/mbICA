#ifndef MBICA_H
#define MBICA_H

#include <armadillo>

//will move to cpp in future
using namespace arma;

namespace mbica {

    const double DEF_EPSILON = 0.0001;
    const int DEF_MAX_ITER = 10000;

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

    class PCA {
    public:
        void operator()(const arma::mat &X, arma::mat &E, arma::vec &D) {
            mat C = cov(X.t());

            eig_sym(D, E, C);

            unsigned stopIx = D.n_elem, startIx=0;
            for(; stopIx > 0 && D[stopIx-1] == 0; --stopIx);
            // Wystarczy sprawdzić od 0 do stopIx, bo w anomalnym przypadku jak wszystkie są 0 to skończymy bez sensu z startIx=D.n_elem i stopIx=0
            for(; startIx < stopIx && D[startIx] == 0; ++startIx);

            if(stopIx != D.n_elem || startIx != 0) {
                D = D.subvec(startIx, stopIx-1);
                // TODO: trzeba sprawdzic, czy eigenvectory w wierszach, czy kolumnach
                // Z tego, co ja rozumiem, powinny być w kolumnach
                E = E.cols(startIx, stopIx-1);
            }

            // TODO: sprawdzic, czy nie potrzebna jest transpozycja E - wszak transponujemy na poczatku X
            // Od transpozycji tutaj zależałoby tylko czy eigenvectory w wyniku będą wierszach, czy kolumnach
        }
    };

    class Whitening {
    public:
        // E - eigenvec from PCA
        // D - vector with eigenvalues from PCA
        void operator()(const arma::mat &E, const arma::vec &D,
                        arma::mat &Wh, arma::mat &dWh) {
            mat sqrtD = sqrt(diagmat(D));
            Wh = inv(sqrtD) * E.t();
            dWh = E * sqrtD;
        }
    };

    mat orth(const mat& A){
        mat U, V;
        vec s;
        svd(U,s,V,A);

        if(A.n_rows > A.n_cols){
            U = U.cols(0, A.n_cols-1);
            s = s.subvec(0, A.n_cols-1);
        }
        double tolerance = std::max(A.n_rows,A.n_cols) * math::eps() * max(s);
        int r;
        for(int i=0; i < s.n_elem; ++i)
            if (s[i] > tolerance)
                ++r;

        return U.cols(0,r-1);
    }

    template<int mu = 1, class UsedFunc = Pow3>
    class FastICA {
    public:
        FastICA(double epsilon = DEF_EPSILON, int maxIterations = DEF_MAX_ITER)
            : whiteningMatrixIsSet_(false) {
            init(epsilon, maxIterations);
        }

        FastICA(mat Wh, mat dWh) {
            init();
            setWhiteningMatrix(Wh, dWh);
        }

        //  tu w jakis posob trzeba umozliwic podanie guess, whitening i dewhitenign matrix
        arma::mat operator()(arma::mat X, int nIC = -1) {            
            //nIC == -1 means the same size it's now.
            if(nIC == -1) {
                nIC = X.n_rows;
            }

            if(!whiteningMatrixIsSet_) {
                mat E;
                vec D;

                PCA()(X, E, D);
                Whitening()(E ,D, Wh_, dWh_);
            }

            // Tu whitening (ktory zawiera w sobie PCA)
            mat A = zeros<mat>(X.n_rows, nIC);
            // B mozemy dac jako zgadniete, np, zeby znalezc wiecej IC
            // B = whiteningMatrix * guess;
            mat B = orth(randu<mat>(X.n_rows, nIC) - 0.5),
                B_old = arma::zeros<mat>(X.n_rows, nIC);
            int i;

            // main loop (easiest way - no stabilization, nor fine-tunung)
            for(i = 0; i < maxIterations_; ++i) {
                B = B * arma::real(matSqrt(arma::inv(B.t() * B)));

                arma::mat minAbsCos = arma::min(arma::abs(arma::diagmat(B.t() * B_old)));
                if( 1 - minAbsCos < epsilon_) {
                    break;
                }
                B_old = B;
                double k = 1.0 / X.n_cols;
                arma::mat Y = X.t() * B;
                B = k * X * UsedFunc::G(Y) - k * arma::repmat(arma::sum(UsedFunc::dG(Y)), X.n_rows, 1) % B;

            }

            if(i == maxIterations_) {
                // return empty A and W
            }
            A = dWh_ * B;
            W = B.t() * Wh_;

            // TODO: returnung std::pair, or maybe by references in parameter list?
            return A, W;
        }

        void setWhiteningMatrix(mat Wh, mat dWh) {
            whiteningMatrixIsSet_ = (dWh != 0 && Wh != 0);
            if(whiteningMatrixIsSet_){
                dWh_ = dWh;
                Wh_ = Wh;
            }
        }

        void setEpsilon(double epsilon) {
            epsilon_ = epsilon;
        }

        void setMaxIterations(int max) {
            maxIterations_ = max;
        }

    private:

        void init(double epsilon = DEF_EPSILON, int maxIterations = DEF_MAX_ITER) {
            epsilon_ = epsilon;
            maxIterations_ = maxIterations;
        }

        bool whiteningMatrixIsSet_;
        mat dWh_;
        mat Wh_;
        double epsilon_;
        int maxIterations_;
    };
}
#endif // MBICA_H
