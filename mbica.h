#ifndef MBICA_H
#define MBICA_H

#include <armadillo>
#include <nonlinearities.h>
#include <icaseparator.h>

//will move to cpp in future
using namespace arma;

namespace mbica {

    const double DEF_EPSILON = 0.0001;
    const int DEF_MAX_ITER = 10000;
    const double DEF_MU = 1.0;

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
        int r=0;
        for(unsigned i=0; i < s.n_elem; ++i)
            if (s[i] > tolerance)
                ++r;

        return U.cols(0,r-1);
    }

    template<class UsedNonl = nonlinearities::Pow<3> >
    class FastICA {
    public:
        FastICA(double epsilon = DEF_EPSILON, int maxIterations = DEF_MAX_ITER, double mu = DEF_MU)
            : whiteningMatrixIsSet_(false) {
            init(epsilon, maxIterations, mu);
        }

        FastICA(mat Wh, mat dWh) {
            init();
            setWhiteningMatrix(Wh, dWh);
        }

        ICASeparator operator()(arma::mat X, int nIC = -1);
        //  tu w jakis posob trzeba umozliwic podanie guess, whitening i dewhitenign matrix

        void setWhiteningMatrix(mat Wh, mat dWh) {
            //whiteningMatrixIsSet_ = (dWh() != 0 && Wh != 0);
            whiteningMatrixIsSet_ = true;
            if(whiteningMatrixIsSet_){
                dWh_ = dWh;
                Wh_ = Wh;
            }
        }

        void unsetWhiteningMatrix() {
            whiteningMatrixIsSet_ = false;
        }

        void setEpsilon(double epsilon) {
            epsilon_ = epsilon;
        }

        void setMaxIterations(int max) {
            maxIterations_ = max;
        }

    private:

        void init(double epsilon = DEF_EPSILON, int maxIterations = DEF_MAX_ITER, double mu = DEF_MU) {
            epsilon_ = epsilon;
            maxIterations_ = maxIterations;
            mu_ = mu;
            stroke_ = 0.0;
            reducedStep_ = false;
        }

        void stabilize(int iteration, arma::mat B, arma::mat B_older ){
            if (stroke_!=0.0)
                mu_=stroke_;
            else{
                double minAbsCos2 = arma::min(arma::min(arma::abs(arma::diagmat(B.t() * B_older))));
                if (1 - minAbsCos2 < epsilon_){
                    stroke_ = mu_;
                    mu_ /= 2;
                }
                else if (!reducedStep_ && (iteration > maxIterations_/2)){
                    reducedStep_=true;
                    mu_ /= 2;
                }
            }
        }

        bool whiteningMatrixIsSet_;
        bool stabilizationEnabled_;
        mat dWh_;
        mat Wh_;
        double epsilon_;
        int maxIterations_;
        double mu_;
        double stroke_;
        bool reducedStep_;
    };

    template<class UsedNonl>
    ICASeparator FastICA<UsedNonl>::operator()(arma::mat X, int nIC) {
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
            B_old = arma::zeros<mat>(X.n_rows, nIC),
            B_older = arma::zeros<mat>(X.n_rows, nIC);
        int i;

        // main loop (with stabilization, but no fine-tunung)
        for(i = 0; i < maxIterations_; ++i) {
            B = B * matSqrt(arma::inv(B.t() * B));

            double minAbsCos = arma::min(arma::min(arma::abs(arma::diagmat(B.t() * B_old))));
            if( 1 - minAbsCos < epsilon_) {
                break;
            }
            else if(stabilizationEnabled_)
                stabilize(i, B, B_older);

            B_older = B_old;
            B_old = B;
            double k = 1.0 / X.n_cols;
            arma::mat Y = X.t() * B;
            UsedNonl nl(Y);
            if(mu_==1.0)
                B = k * X * nl.G() - k * arma::repmat(arma::sum(nl.dG()), X.n_rows, 1) % B;
            else{
                arma::mat Beta = sum(Y % nl.G());
                arma::mat D = diagmat(1/(Beta-sum(nl.dG())));
                B = B + mu_ * B * (Y.t() * nl.G() - diagmat(Beta))* D;
            }


        }

        if(i == maxIterations_) {
            // return empty A and W
        }
        A = dWh_ * B;
        arma::mat W = B.t() * Wh_;

        return ICASeparator(A, W);
    }
}

#endif // MBICA_H
