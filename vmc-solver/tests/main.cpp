#include "UnitTest++/UnitTest++.h"
#include <string.h>
#include <stdio.h>
#include "../system.h"
#include <iostream>
#include <fstream>
#include "../system.h"
#include "../sampler.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/repulsivegaussian.h"
#include "../WaveFunctions/twoelectrons.h"
#include "../WaveFunctions/manyelectrons.h"
#include "../WaveFunctions/manyelectrons_coefficients.h"
#include "../Hamiltonians/squarewell.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../Hamiltonians/harmonicoscillatorrepulsive.h"
#include "../Hamiltonians/harmonicoscillatorelectrons.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../VariationMethods/steepestdescent.h"
#include "../Math/random.h"
#include <mpi.h>
#include <cassert>

using namespace std;
using namespace UnitTest;


SUITE(TestTheTest) {
    TEST(Sanity) {CHECK_EQUAL(1,1);}
}

SUITE(QD) {
    int my_rank = 0;
    int numprocs = 4;
    double totalE, totalKE, totalPE, totalVariance, totalAcceptanceRate, finalMeanDistance;
    double timeStart, timeEnd, totalTime;

    vec L(3);
    L.fill(0.);

    int numberOfDimensions  = 3;
    int numberOfParticles   = 2;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.;
    double alpha            = 0.98456;
    double beta             = 0.40691;
    double gamma            = 2.82843;
    double a                = 0.0043;
    double stepLength       = 0.5;
    double equilibration    = 0.1;
    double dt               = 0.01;
    double aElectrons       = 1.;
    double C                = 1.;
    bool analyticalKinetic  = false;
    bool importanceSampling = true;
    bool repulsion          = true;
    bool quantumDots        = true;
    bool twobodyQD          = false;
    bool Jastrow            = true;
    bool optimizeParameters = false;
    bool saveEnergies       = false;
    bool savePositions      = false;
    bool showProgress       = true;
    bool printToTerminal    = true;
    bool useCoeff 		    = false;

    int numMyCycles = numberOfSteps/numprocs;
    System* system = new System();



    TEST(Initiate) {
        // Initialize MPI parallelization
        MPI_Init(0,0);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        timeStart = MPI_Wtime();

        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, my_rank));

        if (repulsion && !quantumDots) {
        system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, omega, a, gamma, analyticalKinetic));
        system->setWaveFunction             (new RepulsiveGaussian(system, alpha, beta, a));
        }
        if (!repulsion && !quantumDots) {
        system->setHamiltonian              (new HarmonicOscillator(system, omega, analyticalKinetic));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        }
        if (quantumDots) {
            if (twobodyQD) {
                system->setHamiltonian      (new HarmonicOscillatorElectrons(system, omega, analyticalKinetic, repulsion));
                system->setWaveFunction     (new TwoElectrons(system, alpha, beta, omega, aElectrons, C, Jastrow));
            }
            else {
                system->setHamiltonian      (new HarmonicOscillatorElectrons(system, omega, analyticalKinetic, repulsion));
                if (useCoeff) {
                    system->setWaveFunction (new ManyElectronsCoefficients(system, alpha, beta, omega, C, Jastrow));
                }
                else {
                    system->setWaveFunction (new ManyElectrons(system, alpha, beta, omega, C, Jastrow));
                }
            }
        }
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setTimeStep                 (dt);
        system->setMyRank                   (my_rank);
        system->setNumProcs                 (numprocs);

        if (optimizeParameters) {
            system->optimizeParameters          (system, alpha, beta);
        }
        system->setSaveEnergies             (saveEnergies);
        system->setSavePositions            (savePositions);

        system->runMetropolisSteps          (numMyCycles, importanceSampling, showProgress, printToTerminal);

        system->MPI_CleanUp                 (totalE, totalKE, totalPE, totalVariance, totalAcceptanceRate, finalMeanDistance,
                                             timeStart, timeEnd, totalTime, numprocs, numberOfSteps);

        system->mergeOutputFiles            (numprocs);

    }

    TEST(Square) {
        SquareWell* sqwell = new SquareWell(system, 1., 3., omega, false, false);
        /* CHECK(sqwell->computeLocalEnergy(system->getParticles())); */
    }
}



int main() {
    return RunAllTests();
}

