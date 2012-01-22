#define ARMA_USE_LAPACK
#include <armadillo>
#include "mbica.h"

using namespace arma;

int main(){
    mat A;
    A.load("a", raw_ascii);
    mbica::ICASeparator icas = mbica::FastICA<>()(A, -1);

    return 0;
}
