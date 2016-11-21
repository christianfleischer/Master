#include <iostream>
#include <armadillo>
#include <time.h>
#include <cassert>
#include "../system.h"
#include "../Math/factorial.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/doublewell.h"
#include "UnitTest++/UnitTest++.h"
#include "wrapper.h"

using namespace std;
using namespace arma;
using namespace UnitTest;

SUITE(Diagonalization) {
    // Check that factorial function is performing:
    TEST(Factorial) {
        CHECK_EQUAL(87178291200, factorial(14));
    }
    Wrapper* wrapper = new Wrapper();

    TEST(Initialization) {
        int N                       = 1000;
        double posMin               = -10;
        double posMax               = 10;
        double omega_r              = 0.5;                                         // =m*w/hbar Just a constant to keep the results correct, while we figure out the omega conundrum.
        int nMax 					= 3;
        int nPrimeMax               = 2;
        int numberOfDimensions      = 2;

        wrapper->m_N = N;
        wrapper->m_posMin = posMin;
        wrapper->m_posMax = posMax;
        wrapper->m_omega_r = omega_r;
        wrapper->m_nMax = nMax;
        wrapper->m_nPrimeMax = nPrimeMax;
        wrapper->m_numberOfDimensions = numberOfDimensions;

        bool harmOscPotential       = true;

        wrapper->m_harmOscPotential = harmOscPotential;

        vec L(numberOfDimensions);
        L.fill(0.);
        L(0) = 0.;

        wrapper->setL(L);

        int numberOfEigstates;


        if (wrapper->m_numberOfDimensions == 2) {
            //numberOfEigstates = int(0.5*(nMax+1)*(nMax+2));
            numberOfEigstates = int(0.5*(wrapper->m_nMax)*(wrapper->m_nMax+1));
        }

        else if (wrapper->m_numberOfDimensions == 3) {
            //numberOfEigstates = int((nMax+1)*(nMax+2)*(nMax+3)/6.);
            numberOfEigstates = int((wrapper->m_nMax)*(wrapper->m_nMax+1)*(wrapper->m_nMax+2)/6.);
        }
        else { numberOfEigstates = wrapper->m_nMax; }

        assert(wrapper->m_nPrimeMax <= numberOfEigstates);

        wrapper->m_numberOfEigstates = numberOfEigstates;

    }

    TEST(McTestFace) {



        //Set up the vector x and the matrix A:
        double h                    = (wrapper->m_posMax-wrapper->m_posMin)/wrapper->m_N;
        mat r = zeros(wrapper->m_N+1,wrapper->m_numberOfDimensions);
        vec rAbs = zeros(wrapper->m_N+1);

        for (int d = 0; d < wrapper->m_numberOfDimensions; d++) {
            r.col(d) = wrapper->m_posMin + linspace(0, wrapper->m_N, wrapper->m_N+1)*h;
            rAbs += r.col(d)%r.col(d);
        }

        rAbs = sqrt(rAbs);


        cube diagMat(wrapper->m_N-1, wrapper->m_N-1, wrapper->m_numberOfDimensions);
        cube eigvecs(wrapper->m_N-1, wrapper->m_N-1, wrapper->m_numberOfDimensions);
        mat  eigvals(wrapper->m_N-1, wrapper->m_numberOfDimensions);

        mat saveEigenvector         = ones(wrapper->m_N-1, wrapper->m_numberOfEigstates);
        cube saveSepEigenvector		= zeros(wrapper->m_N-1, wrapper->m_numberOfEigstates, wrapper->m_numberOfDimensions);
        mat SavePositionvector      = zeros(wrapper->m_N-1, wrapper->m_numberOfDimensions+1);

        for (int d = 0; d < wrapper->m_numberOfDimensions; d++) {
            SavePositionvector.col(d)   = r.col(d).subvec(1, wrapper->m_N-1);    //Saves the y vector for output.
        }

        SavePositionvector.col(wrapper->m_numberOfDimensions) = rAbs.subvec(1, wrapper->m_N-1);    //Saves the r vector for output.

        //Init system
        System* system = new System(wrapper->m_omega_r, wrapper->m_numberOfDimensions, h, wrapper->m_N);

        system->setWaveFunction(new DoubleWell(system, wrapper->m_omega_r));

        system->diagonalizeMatrix(r, wrapper->m_L, wrapper->m_N, diagMat);
        system->findEigenstate(eigvals, eigvecs, diagMat,
                               saveEigenvector, saveSepEigenvector,
                               wrapper->m_numberOfEigstates, wrapper->m_nMax);

        // Check to see if eigenvalues correspond to the first values for the given potential.
        if (wrapper->m_harmOscPotential) {
            CHECK_CLOSE(0.5, eigvals(0,0), 0.0001);
        }

        cout << endl << "eigvals, Armadillo:" << endl;
        int displayVals = 15;
        for (int i = 0; i < displayVals; ++i) {
            for (int d = 0; d < wrapper->m_numberOfDimensions; d++) {
                cout << i+1 << ": E"<< d <<": " << eigvals.col(d)(i);
            }
            cout << endl;
        }
        cout << endl;
        cout << "Computation time (sec):" << endl;
        cout << "Aramadillo: " << system->getComputationTime() << endl;
    }
}

int main() {
    return RunAllTests();
}
