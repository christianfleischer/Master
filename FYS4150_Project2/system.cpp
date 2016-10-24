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

void System::findEigenstate(mat &eigvals, cube eigvecs, cube diagMat,
                            mat &saveEigenvector,
                            cube &saveSepEigenvector,
                            int numberOfEigstates) {
    clock_t start1, finish1;
    start1 = clock();

    // Finding eigenvalues and eigenvectors using armadillo:
    for (int d = 0; d < m_numberOfDimensions; d++) {
        vec eigvalsTemp = eigvals.col(d);
        eig_sym(eigvalsTemp, eigvecs.slice(d),  diagMat.slice(d));
        eigvals.col(d) = eigvalsTemp;
    }

    cube eigVecsTemp = zeros(m_N-1, numberOfEigstates, m_numberOfDimensions);
    for (int d = 0; d < m_numberOfDimensions; d++) {    // !!!!!!!
        eigVecsTemp.slice(d) = eigvecs.slice(d).submat(0,0,m_N-2,numberOfEigstates-1);
        saveEigenvector %= eigVecsTemp.slice(d);
        saveSepEigenvector.slice(d) = eigVecsTemp.slice(d);
    }

    m_psi = saveSepEigenvector;

    finish1 = clock();
    m_computationTime = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    return;
}

void System::findCoefficients(int nMax, int nPrimeMax, vec x, mat &C, int currentDim){
    cout << "Finding coefficients for dimension " << currentDim+1 << " of " <<  m_numberOfDimensions << endl;
    cout.flush();
    std::string upLine = "\033[F";
    for	(int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
        cout << "nPrime = " << nPrime << " of " << nPrimeMax-1 << endl;
        for (int nx = 0; nx < nMax; nx++) {
            cout << "[" << int(double(nx)/nMax * 100.0) << " %]\r";
            cout.flush();
            double innerprod = 0;
            for (int i = 0; i < m_N-1; i++) {
                innerprod += m_psi.slice(currentDim).col(nPrime)(i)*m_waveFunction->harmonicOscillatorBasis(x, nx)(i);
            }
            C(nx, nPrime) = innerprod;
        }
        cout << upLine;
        //nPrime++;   //Only need even nPrimes due to double well (degeneracy = 2).
    }
    cout << upLine;
    C *= m_h;
}

mat System::findSuperPos(mat r, int nMax, int nPrimeMax, cube &supPosSep, mat &saveC) {
    mat rCut = zeros(m_N-1, m_numberOfDimensions);

    for (int d=0; d < m_numberOfDimensions; d++) {
        rCut.col(d) = r.col(d).subvec(1,m_N-1);
    }

    //rCut.col(1) = zeros(m_N-1);
    //rCut.col(2) = zeros(m_N-1);

    cube C = zeros(nMax, nPrimeMax, m_numberOfDimensions);
    mat supPos = ones(m_N-1, nPrimeMax);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        findCoefficients(nMax, nPrimeMax, rCut.col(d), C.slice(d), d);
        saveC %= C.slice(d);
    }

//    for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
//        for (int n = 0; n < nMax; n++) {
//            vec plusTerm = ones(m_N-1);
//            for (int d = 0; d < m_numberOfDimensions; d++) {
//                plusTerm %= C(n, nPrime, d)*m_waveFunction->harmonicOscillatorBasis(rCut.col(d), n);
//            }
//            supPos.col(nPrime) += plusTerm;
//        }
//        nPrime++;   //Only need even nPrimes due to double well (degeneracy = 2).
//    }

    for (int nPrime = 0; nPrime < nPrimeMax; nPrime++) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            vec plusTerm = zeros(m_N-1);
            for (int n = 0; n < nMax; n++) {
                plusTerm += C(n, nPrime, d)*m_waveFunction->harmonicOscillatorBasis(rCut.col(d), n);
            }
            supPos.col(nPrime) %= plusTerm;
            supPosSep.slice(d).col(nPrime) = plusTerm;
        }
        //nPrime++;   //Only need even nPrimes due to double well (degeneracy = 2).
    }

    return supPos;
}


void System::setWaveFunction(WaveFunction *waveFunction) { m_waveFunction = waveFunction; }
void System::setN(double N) {m_N = N; }
void System::setStepLength(double h) { m_h = h; }
void System::setNumberOfDimensions(int numberOfDimensions) { m_numberOfDimensions = numberOfDimensions; }
void System::setNumberOfParticles(int numberOfParticles) { m_numberOfParticles = numberOfParticles; }
void System::setComputationTime(double computationTime) { m_computationTime = computationTime; }
