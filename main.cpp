#include <armadillo>
#include "mbica.h"

using namespace arma;
using namespace mbica;

int main(){
    mat A, guess;
    A.load("a", raw_ascii);
    guess.load("g", raw_ascii);
    ICASeparator icas = FastICA<>(guess)(A, -1);

    return 0;
}
