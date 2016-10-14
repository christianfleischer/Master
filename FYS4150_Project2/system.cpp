#include "system.h"
#include "WaveFunctions/wavefunction.h"

System::System(double omega, int numberOfDimensions, double h, int N) {
    m_omega = omega;
    setNumberOfDimensions(numberOfDimensions);
    setN(N);
    setStepLength(h);
}

void System::diagonalizeMatrix(mat r, vec L, int N, cube &diagMat) {
    double Constant = 1./(m_h*m_h);
    mat V(N+1, m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {
        V.col(d) = m_waveFunction->potential(r.col(d), L(d));
        diagMat.slice(d).diag(0)  =  2*Constant + V.col(d).subvec(1,N-1);     //Set d_i elements in A
        diagMat.slice(d).diag(1)  = -Constant*ones(N-2);               //Set e_i elements in A
        diagMat.slice(d).diag(-1) = diagMat.slice(d).diag(1);                         //Set e_i elements in A
    }

    return;
}

void System::findEigenstate(mat &eigvals, cube eigvecs, cube diagMat, int numberOfEigstates, mat &saveEigenvector) {
    clock_t start1, finish1;
    start1 = clock();

    // Finding eigenvalues and eigenvectors using armadillo:
    for (int d = 0; d < m_numberOfDimensions; d++) {
        vec eigvalsTemp = eigvals.col(d);
        eig_sym(eigvalsTemp, eigvecs.slice(d),  diagMat.slice(d));
        eigvals.col(d) = eigvalsTemp;
    }

    for (int i = 0; i < numberOfEigstates; i++) {
        for (int d = 0; d < 1; d++) {
            saveEigenvector.col(i) %= eigvecs.slice(d).col(i);
        }
    }

    //m_waveFunction.setWaveFunction(saveEigenvector);

    finish1 = clock();
    m_computationTime = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    return;
}

void System::findCoefficients(mat r, vec qNumbers, mat &C){
    //C = waveFunction.getWaveFunction().col(0)%waveFunction.harmonicOscillatorBasis(r, qNumbers);
}


void System::setWaveFunction(WaveFunction *waveFunction) { m_waveFunction = waveFunction; }
void System::setN(double N) {m_N = N; }
void System::setStepLength(double h) { m_h = h; }
void System::setNumberOfDimensions(int numberOfDimensions) { m_numberOfDimensions = numberOfDimensions; }
void System::setNumberOfParticles(int numberOfParticles) { m_numberOfParticles = numberOfParticles; }
void System::setComputationTime(double computationTime) { m_computationTime = computationTime; }
