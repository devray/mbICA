#include <armadillo>
#include <cstdlib>
#include <iostream>
#include "mbica.h"

using namespace std;
using namespace arma;
using namespace mbica;

int main(int argc, char *argv[]){
    if(argc < 2) {
        return -1;
    }
    double mu = 1.0;
    if(argc > 2) {
        mu = atof(argv[2]);
    }

    cout << "Using mu = " << mu << endl;
    mat A;
    A.load(argv[1], raw_ascii);
    A = A.t();
    ICASeparator icas = FastICA<nonlinearities::Pow<3>, WithStabilization>(_mu = mu)(A,15);
    mat(icas(A).t()).save("result.txt", raw_ascii);

    return 0;
}
