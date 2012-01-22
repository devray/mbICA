#define ARMA_USE_LAPACK
#include <armadillo>
#include "mbica.h"

using namespace arma;

int main(){
    mat A, guess;
    A.load("a", raw_ascii);
    guess.load("g", raw_ascii);
    mbica::ICASeparator icas = mbica::FastICA<>(guess)(A, -1);

    return 0;
}
