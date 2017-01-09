#include <iostream>
#include <armadillo>
#include <time.h>
#include <cassert>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/doublewell.h"
#include "WaveFunctions/finitewell.h"
#include "WaveFunctions/squarewell.h"

using namespace std;
using namespace arma;

int main() {

    int N                       = 1000;
    double posMin               = -10;
    double posMax               = 10;
    double omega_r              = 0.5;                                         // =m*w/hbar Just a constant to keep the results correct, while we figure out the omega conundrum.
    double V0                   = 1.;
    int nMax 					= 4;
    int nPrimeMax               = 4;
    int numberOfDimensions      = 3;
    double distanceToWall       = 3.;

    vec L(3);
    L.fill(0.);
    L(0) = 5.;

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
    cube supPosSep				= zeros(N-1, nPrimeMax, numberOfDimensions);
    cube saveC = ones(nMax, nPrimeMax, numberOfDimensions);

    for (int d = 0; d < numberOfDimensions; d++) {
        SavePositionvector.col(d)   = r.col(d).subvec(1, N-1);    //Saves the y vector for output.
    }

    SavePositionvector.col(numberOfDimensions) = rAbs.subvec(1, N-1);    //Saves the r vector for output.

    vec SaveConstants           = {omega_r, double(numberOfDimensions), L(0), L(1), L(2), double(N), double(numberOfEigstates), h};

    //Init system
    System* system = new System(omega_r, numberOfDimensions, h, N);

    system->setWaveFunction(new DoubleWell(system, omega_r));
    //system->setWaveFunction(new FiniteWell(system, omega_r, distanceToWall));
    //system->setWaveFunction(new SquareWell(system, omega_r, V0, distanceToWall));

    system->diagonalizeMatrix(r, L, N, diagMat);
    system->findEigenstate(eigvals, eigvecs, diagMat,
                           saveEigenvector, saveSepEigenvector,
                           numberOfEigstates, nMax);
    mat supPos = system->findSuperPos(r, nMax, nPrimeMax, supPosSep, saveC);

    SaveConstants.save("../diagonalization/PlotAndData/Constants.dat", raw_ascii);
    SavePositionvector.save("../diagonalization/PlotAndData/Positionvectors.dat", raw_ascii);
    saveEigenvector.save("../diagonalization/PlotAndData/Eigenvectors.dat", raw_ascii);
    saveSepEigenvector.slice(0).save("../diagonalization/PlotAndData/SeparateEigenvectorsX.dat", raw_ascii);
    if (numberOfDimensions>1) saveSepEigenvector.slice(1).save("../diagonalization/PlotAndData/SeparateEigenvectorsY.dat", raw_ascii);
    supPos.save("../diagonalization/PlotAndData/Superpositions.dat", raw_ascii);
    supPosSep.slice(0).save("../diagonalization/PlotAndData/SeparateSuperpositionsX.dat", raw_ascii);
    if (numberOfDimensions>1) supPosSep.slice(1).save("../diagonalization/PlotAndData/SeparateSuperpositionsY.dat", raw_ascii);
    saveC.save("../diagonalization/PlotAndData/Coefficients.dat", arma_ascii);
    //saveEigenvector.print();
    eigvals.save("../diagonalization/PlotAndData/Eigenvalues.dat", arma_ascii);

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

    return 0;
}
