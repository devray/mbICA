#ifndef ICASEPARATOR_H
#define ICASEPARATOR_H

#include <armadillo>

namespace mbica {
class ICASeparator
{
public:
    ICASeparator(arma::mat A, arma::mat W)
        : A_(A), W_(W) {}

    arma::mat operator()(arma::mat X);

private:
    arma::mat A_;
    arma::mat W_;
};
}

#endif // ICASEPARATOR_H
