#include <iostream>
#include "quantumdots.h"

using namespace std;

int main()
{
    int numberOfShells = 4;
    int numberOfParticles = 6;
    int numberOfDimensions = 2;
    double omega = 1.;

    QuantumDots quantumDots(numberOfShells, numberOfParticles, numberOfDimensions, omega);
    quantumDots.runHartreeFock();
}

