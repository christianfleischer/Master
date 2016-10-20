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

    m_psi = saveEigenvector;

    finish1 = clock();
    m_computationTime = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    return;
}

void System::findCoefficients(mat r, vec qNumbers, vec &C){
    int nMax = 18;
    for (int nx = 0; nx < nMax; nx++) {
        double innerprod = 0;
        vec promp = {nx,0,0};
        for (int i = 0; i < m_N-1; i++) {
            innerprod += m_h*m_psi.col(0)(i)*m_waveFunction->harmonicOscillatorBasis(r, promp)(i);
        }
        C(nx) = innerprod;
    }
    //C = m_psi.col(0)%m_waveFunction->harmonicOscillatorBasis(r, qNumbers);
}

vec System::findSuperPos(mat r, int nMax) {
    mat rCut = zeros(m_N-1, m_numberOfDimensions);

    for (int d=0; d < m_numberOfDimensions; d++) {
        rCut.col(d) = r.col(d).subvec(1,m_N-1);
    }

    //rCut.col(1) = zeros(m_N-1);
    //rCut.col(2) = zeros(m_N-1);

    vec C = zeros(nMax);
    vec supPos = zeros(m_N-1);
  //  cout << "Starting supPos calculation." << endl;
  //  for (int nx = 0; nx < nMax; nx++) {
  //      cout << double(nx)/nMax*100 << endl;
  //      for (int ny = 0; ny < nMax; ny++) {
  //          for (int nz = 0; nz < nMax; nz++) {
  //              vec qNumbers = {double(nx), double(ny), double(nz)};
  //              findCoefficients(rCut, qNumbers, C);
  //              vec plusTerm = C%m_waveFunction->harmonicOscillatorBasis(rCut, qNumbers);
  //              supPos += plusTerm;
  //          }
  //      }
  //  }


    findCoefficients(rCut, {0,0,0}, C);
    for (int nx = 0; nx < nMax; nx++) {
        vec qNumbers = {double(nx),0,0};
        vec plusTerm = C(nx)*m_waveFunction->harmonicOscillatorBasis(rCut, qNumbers);
        supPos += plusTerm;
    }

    return supPos;
}


void System::setWaveFunction(WaveFunction *waveFunction) { m_waveFunction = waveFunction; }
void System::setN(double N) {m_N = N; }
void System::setStepLength(double h) { m_h = h; }
void System::setNumberOfDimensions(int numberOfDimensions) { m_numberOfDimensions = numberOfDimensions; }
void System::setNumberOfParticles(int numberOfParticles) { m_numberOfParticles = numberOfParticles; }
void System::setComputationTime(double computationTime) { m_computationTime = computationTime; }
