#include <armadillo>
#include "mbica.h"

using namespace arma;
using namespace mbica;

int main(int argc, char *argv[]){
    if(argc < 2) {
        return -1;
    }
    mat A;
    A.load(argv[1], raw_ascii);
    A = A.t();
    ICASeparator icas = FastICA<>()(A);//_guessMatrix = guess)(A);
    mat(icas(A).t()).save("result.txt", raw_ascii);

    return 0;
}
