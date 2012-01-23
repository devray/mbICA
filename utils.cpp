#include "utils.h"

using namespace arma;

namespace mbica {

arma::mat matSqrt(const arma:: mat &X) {
    vec D;
    mat E;
    eig_sym(D, E, X);

    // TODO: Sprawdzic, czy znormalizowane eigenvectory

    return E * diagmat(sqrt(D)) * E.t();
}

void PCA::operator()(const arma::mat &X, arma::mat &E, arma::vec &D) {
    mat C = cov(X.t());

    eig_sym(D, E, C);

    unsigned stopIx = D.n_elem, startIx=0;
    for(; stopIx > 0 && D[stopIx-1] < math::eps(); --stopIx);
    // Wystarczy sprawdzić od 0 do stopIx, bo w anomalnym przypadku jak wszystkie są 0 to skończymy bez sensu z startIx=D.n_elem i stopIx=0
    for(; startIx < stopIx && D[startIx] < math::eps(); ++startIx);

    if(stopIx != D.n_elem || startIx != 0) {
        D = D.subvec(startIx, stopIx-1);
        // TODO: trzeba sprawdzic, czy eigenvectory w wierszach, czy kolumnach
        // Z tego, co ja rozumiem, powinny być w kolumnach
        E = E.cols(startIx, stopIx-1);
    }

    // TODO: sprawdzic, czy nie potrzebna jest transpozycja E - wszak transponujemy na poczatku X
    // Od transpozycji tutaj zależałoby tylko czy eigenvectory w wyniku będą wierszach, czy kolumnach
}
void Whitening::operator()(const arma::mat &E, const arma::vec &D,
                           arma::mat &Wh, arma::mat &dWh) {
    mat sqrtD = diagmat(sqrt(D));
    Wh = inv(sqrtD) * E.t();
    dWh = E * sqrtD;
}

mat orth(const mat &A){
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

mat remmean(const mat &X) {
    return X - arma::repmat(arma::mean(X, 1), 1, X.n_cols);
}
}
