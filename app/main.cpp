#include <armadillo>
#include <cstdlib>
#include <iostream>
#include <string>
#include "mbica.h"

using namespace std;
using namespace arma;
using namespace mbica;

int main(int argc, char *argv[]){
    if(argc < 2)
        return -1;

    double mu = 1.0;
    if(argc > 2)
        mu = atof(argv[2]);

    // Setting rand seed
    srand(time(NULL));

    cout << "Using mu = " << mu << endl;
    mat A;
    A.load(argv[1], raw_ascii);

    // We have to remove two central signals as they are not EEG signals
    A = join_rows(A.cols(0, 6), A.cols(9, 15));

    // Signals are in columns so we have to transpose to get signals in rows
    A = A.t();
    ICASeparator icas = FastICA<nonlinearities::Gauss<>, WithStabilization>()(A, -1, mu);

    // Saving results
    mat(icas(A).t()).save("result_" + string(argv[1]), raw_ascii);

    return 0;
}
