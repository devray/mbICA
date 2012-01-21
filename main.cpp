#include <armadillo>
#include "mbica.h"

using namespace arma;

int main(){
    mat A;
    mbica::FastICA<>()(A);

    return 0;
}
