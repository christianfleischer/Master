#include <iostream>
#include <armadillo>
#include <time.h>
#include <cassert>
#include "../system.h"
#include "../Math/factorial.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/doublewell.h"
#include "UnitTest++/UnitTest++.h"

using namespace std;
using namespace arma;
using namespace UnitTest;

void initSystem() {

}

TEST(Diagonalization) {
    int N                       = 1000;
    double posMin               = -10;
    double posMax               = 10;
    double omega_r              = 0.5;                                         // =m*w/hbar Just a constant to keep the results correct, while we figure out the omega conundrum.
    int nMax 					= 3;
    int nPrimeMax               = 2;
    int numberOfDimensions      = 2;

    bool harmOscPotential       = true;

    vec L(3);
    L.fill(0.);
    L(0) = 0.;

    int numberOfEigstates;

    if (numberOfDimensions == 2) {
        //numberOfEigstates = int(0.5*(nMax+1)*(nMax+2));
        numberOfEigstates = int(0.5*(nMax)*(nMax+1));
    }

    else if (numberOfDimensions == 3) {
        //numberOfEigstates = int((nMax+1)*(nMax+2)*(nMax+3)/6.);
        numberOfEigstates = int((nMax)*(nMax+1)*(nMax+2)/6.);
    }
    else { numberOfEigstates = nMax; }

    assert(nPrimeMax <= numberOfEigstates);

    // Check that factorial function is performing:
    CHECK_EQUAL(87178291200, factorial(14));

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
    cube saveSepEigenvector		= zeros(N-1, numberOfEigstates, numberOfDimensions);
    mat SavePositionvector      = zeros(N-1, numberOfDimensions+1);

    for (int d = 0; d < numberOfDimensions; d++) {
        SavePositionvector.col(d)   = r.col(d).subvec(1, N-1);    //Saves the y vector for output.
    }

    SavePositionvector.col(numberOfDimensions) = rAbs.subvec(1, N-1);    //Saves the r vector for output.

    //Init system
    System* system = new System(omega_r, numberOfDimensions, h, N);

    system->setWaveFunction(new DoubleWell(system, omega_r));

    system->diagonalizeMatrix(r, L, N, diagMat);
    system->findEigenstate(eigvals, eigvecs, diagMat,
                           saveEigenvector, saveSepEigenvector,
                           numberOfEigstates, nMax);

    // Check to see if eigenvalues correspond to the first values for the given potential.
    if (harmOscPotential) {
        CHECK_CLOSE(0.5, eigvals(0,0), 0.0001);
    }

    cout << endl << "eigvals, Armadillo:" << endl;
    int displayVals = 15;
    for (int i = 0; i < displayVals; ++i) {
        for (int d = 0; d < numberOfDimensions; d++) {
            cout << i+1 << ": E"<< d <<": " << eigvals.col(d)(i);
        }
        cout << endl;
    }
    cout << endl;
    cout << "Computation time (sec):" << endl;
    cout << "Aramadillo: " << system->getComputationTime() << endl;


}

int main() {
    return RunAllTests();
}
