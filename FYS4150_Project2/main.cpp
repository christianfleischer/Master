#include <iostream>
#include <armadillo>
#include <time.h>
#include <potential.cpp>
#include <initialize_electrons.cpp>


using namespace std;
using namespace arma;

int main()
{
    int N = 1000;
    double posMin = -10;
    double posMax = 10;
    double omega_r = 0.5;                                         // =m*w/hbar Just a constant to keep the results correct, while we figure out the omega conundrum.

    vec L(3);
    L.fill(5.);
    mat X = zeros(N-1,N-1);
    mat Y = zeros(N-1,N-1);
    mat Z = zeros(N-1,N-1);

    int numEigvectors = 100;

    int nDim = 3;

    //Set up the vector x and the matrix A:
    double h = (posMax-posMin)/N;
    vec x = posMin + linspace(0, N, N+1)*h;
    vec y = posMin + linspace(0, N, N+1)*h;
    vec z = posMin + linspace(0, N, N+1)*h;
    //vec rho = linspace(0, N, N+1)*h;

    vec rho = sqrt(x%x + y%y + z%z);

    vec Vx, Vy, Vz;

    Potential(Vx, x, L(0), omega_r);
    Potential(Vy, y, L(1), omega_r);
    Potential(Vz, z, L(2), omega_r);


    mat SaveEigenvector = zeros(N-1, numEigvectors);

    mat SavePositionvector = zeros(N-1, nDim+1);
    SavePositionvector.col(0) = x.subvec(1, N-1);    //Saves the x vector for output.
    SavePositionvector.col(1) = y.subvec(1, N-1);    //Saves the y vector for output.
    SavePositionvector.col(2) = z.subvec(1, N-1);    //Saves the z vector for output.
    SavePositionvector.col(3) = rho.subvec(1, N-1);    //Saves the r vector for output.

    vec SaveConstants = {omega_r, L(0), L(1), L(2), double(N), double(numEigvectors)};

    //Set initial condtions and set up matrix:
    InitializeOneElectron(N, X, Vx, h);
    InitializeOneElectron(N, Y, Vy, h);
    InitializeOneElectron(N, Z, Vz, h);

    clock_t start1, finish1;
    start1 = clock();

    //Finding eigenvalues and eigenvectors using armadillo:
    vec eigvalsX;
    vec eigvalsY;
    vec eigvalsZ;

    mat eigvecsX;
    mat eigvecsY;
    mat eigvecsZ;

    eig_sym(eigvalsX, eigvecsX,  X);
    eig_sym(eigvalsY, eigvecsY,  Y);
    eig_sym(eigvalsZ, eigvecsZ,  Z);

    finish1 = clock();
    double ComputationTimeArma = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    for (int i = 0; i < numEigvectors; ++i) {
        SaveEigenvector.col(i) = eigvecsX.col(i)%eigvecsY.col(i)%eigvecsZ.col(i);
    }

    SaveConstants.save("/home/alexanfl/master/FYS4150_Project2/PlotAndData/omega5constants.dat", raw_ascii);
    SavePositionvector.save("/home/alexanfl/master/FYS4150_Project2/PlotAndData/omega5position.dat", raw_ascii);
    SaveEigenvector.save("/home/alexanfl/master/FYS4150_Project2/PlotAndData/omega5norepulsion.dat", raw_ascii);

    cout << "eigvals, Armadillo:" << endl;
    int displayVals = 15;
    for (int i = 0; i < displayVals; ++i) {
        cout << i+1 << ": Ex: " << eigvalsX(i) << " Ey: " << eigvalsY(i) << " Ez: " << eigvalsZ(i) << endl;
    }
    cout << endl;
    cout << "Computation time (sec):" << endl;
    cout << "Aramadillo: " << ComputationTimeArma << endl;

    return 0;
}
