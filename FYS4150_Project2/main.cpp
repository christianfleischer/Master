#include <iostream>
#include <armadillo>
#include <time.h>
#include "system.h"

using namespace std;
using namespace arma;

int main() {
    vec L(3);
    L.fill(5.);

    int N                       = 100;
    double posMin               = -10;
    double posMax               = 10;
    double omega_r              = 0.5;                                         // =m*w/hbar Just a constant to keep the results correct, while we figure out the omega conundrum.

    int numberOfEigstates       = 50;
    int numberOfDimensions      = 3;

    //Set up the vector x and the matrix A:
    double h                    = (posMax-posMin)/N;
    mat r = zeros(N+1,numberOfDimensions);
    vec rAbs = zeros(N+1);

    for (int d = 0; d < numberOfDimensions; d++) {
        r.col(d) = posMin + linspace(0, N, N+1)*h;
        rAbs += r.col(d)%r.col(d);
    }
    rAbs = sqrt(rAbs);


    cube diagMat(N-1, N-1, numberOfDimensions);
    cube eigvecs(N-1, N-1, numberOfDimensions);
    mat  eigvals(N-1, numberOfDimensions);

    mat saveEigenvector         = ones(N-1, numberOfEigstates);

    mat SavePositionvector      = zeros(N-1, numberOfDimensions+1);

    for (int d = 0; d < numberOfDimensions; d++) {
        SavePositionvector.col(d)   = r.col(d).subvec(1, N-1);    //Saves the y vector for output.

    }

    SavePositionvector.col(numberOfDimensions) = rAbs.subvec(1, N-1);    //Saves the r vector for output.

    vec SaveConstants           = {omega_r, L(0), L(1), L(2), double(N), double(numberOfEigstates)};

    //Init system
    System* system = new System(omega_r, numberOfDimensions, h);
    system->diagonalizeMatrix(r, L, N, diagMat);
    system->findEigenstate(eigvals, eigvecs, diagMat, numberOfEigstates, saveEigenvector);

    SaveConstants.save("../FYS4150_Project2/PlotAndData/omega5constants.dat", raw_ascii);
    SavePositionvector.save("../FYS4150_Project2/PlotAndData/omega5position.dat", raw_ascii);
    saveEigenvector.save("../FYS4150_Project2/PlotAndData/omega5norepulsion.dat", raw_ascii);
    //saveEigenvector.print();

    cout << "eigvals, Armadillo:" << endl;
    int displayVals = 15;
    for (int i = 0; i < displayVals; ++i) {
        cout << i+1 << ": Ex: " << eigvals.col(0)(i) << " Ey: " << eigvals.col(1)(i) << " Ez: " << eigvals.col(2)(i) << endl;
    }
    cout << endl;
    cout << "Computation time (sec):" << endl;
    cout << "Aramadillo: " << system->getComputationTime() << endl;

    return 0;
}
