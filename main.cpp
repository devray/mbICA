#define ARMA_USE_LAPACK
#include <armadillo>
#include "mbica.h"

using namespace arma;

int main(){
    mat A=randu<mat>(3,6);
    mbica::ICASeparator icas = mbica::FastICA<>()(A, -1);

    return 0;
}
