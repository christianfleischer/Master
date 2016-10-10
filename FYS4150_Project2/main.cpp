#include <iostream>
#include <armadillo>
#include <time.h>
#include <initialize_electrons.cpp>
#include <jacobi_rotations.cpp>

using namespace std;
using namespace arma;

int main()
{
    int N = 2000;
    double posMin = -10;
    double posMax = 10;
    double omega = 1;                                         // =m*w/hbar Just a constant to keep the results correct, while we figure out the omega conundrum.

    vec L(3);
    L.fill(4.);

    mat A = zeros(N-1,N-1);

    int numEigvectors = 200;

    int nDim = 2;

    //Set up the vector x and the matrix A:
    double h = (posMax-posMin)/N;
    vec x = posMin + linspace(0, N, N+1)*h;
    vec y = posMin + linspace(0, N, N+1)*h;
    vec V = 0.5*omega*(x%x + y%y
                       - 2*abs(x)*L(0) - 2*abs(y)*L(1)
                       + L(0)*L(0) + L(1)*L(1));

    mat SaveEigenvector = zeros(N-1, numEigvectors);
    mat SavePositionvector = zeros(N-1, nDim);
    SavePositionvector.col(0) = x.subvec(1, N-1);    //Saves the x vector for output.
    SavePositionvector.col(1) = y.subvec(1, N-1);    //Saves the y vector for output.

    vec SaveConstants(5);
    SaveConstants(0) = omega;
    SaveConstants(1) = L(0);
    SaveConstants(2) = L(1);
    SaveConstants(3) = N;
    SaveConstants(4) = numEigvectors;


    //Set initial condtions and set up matrix:
    InitializeOneElectron(N, A, omega, V, h);
    //InitializeTwoElectrons(N, A, xMin, xMax, SaveEigenvector, omega_r);

    clock_t start1, finish1;
    start1 = clock();

    //Finding eigenvalues and eigenvectors using armadillo:
    vec ArmadilloEigenvalues;
    mat Eigenvectors;
    eig_sym(ArmadilloEigenvalues, Eigenvectors, A);

    finish1 = clock();
    double ComputationTimeArma = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    //Finding eigenvalues using Jacobi rotations:
    int NumberOfRotations = 0;      //Counter for similarity transformations.
    //double Tolerance = 1E-8;        //Tolerance for when the matrix should be considered as diagonal.
    mat absA = abs(A);              //Creating a matrix for finding the max value element in A.
    absA.diag(0) = zeros(N-1);      //Remove the diagonal elements since we're interested in the max off-diagonal element.

    //uword RowIndexMax;      //Row index for the maximum element.
    //uword ColIndexMax;      //Column index for the maximum element.
    //double MaxValue = absA.max(RowIndexMax, ColIndexMax);   //Value of the max element.

    clock_t start2, finish2;
    start2 = clock();

    //Diagonalize A with Jacobi rotations:
        //JacobiRotations(MaxValue, A, absA, NumberOfRotations, Tolerance, RowIndexMax, ColIndexMax, N);
    //Eigenvalues listed on the diagonal of A sorted from smallest to largest:
    vec JacobiEigenvalues = sort(A.diag(0));

    finish2 = clock();
    double ComputationTimeJacobi = ((finish2-start2)/(double) CLOCKS_PER_SEC);



    for (int i = 0; i < numEigvectors; ++i) {
        SaveEigenvector.col(i) = Eigenvectors.col(i);
    }

    SaveConstants.save("/home/alexanfl/master/FYS4150_Project2/PlotAndData/omega5constants.dat", raw_ascii);
    SavePositionvector.save("/home/alexanfl/master/FYS4150_Project2/PlotAndData/omega5position.dat", raw_ascii);
    SaveEigenvector.save("/home/alexanfl/master/FYS4150_Project2/PlotAndData/omega5norepulsion.dat", raw_ascii);

    cout << "Eigenvalues, Jacobi:" << endl;
    cout << "1:     " << JacobiEigenvalues(0) << endl;
    cout << "2:     " << JacobiEigenvalues(1) << endl;
    cout << "3:     " << JacobiEigenvalues(2) << endl;
    cout << "Eigenvalues, Armadillo:" << endl;

    int displayVals = 15;
    for (int i = 0; i < displayVals; ++i) {
        cout << i+1 << ":     " << ArmadilloEigenvalues(i) << endl;
    }
    cout << endl;
    cout << "Number of similarity transformations: " << NumberOfRotations << endl;
    cout << "Computation time (sec):" << endl;
    cout << "Jacobi: " << ComputationTimeJacobi << "    " << "Aramadillo: " << ComputationTimeArma << endl;

    return 0;
}
