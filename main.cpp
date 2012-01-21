#define ARMA_USE_LAPACK
#include <armadillo>
#include "mbica.h"

using namespace arma;

int main(){
    mat A;
    mbica::FastICA<>()(A, -1);

    return 0;
}
