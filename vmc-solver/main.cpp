#include <iostream>
#include <unittest++/UnitTest++.h>
#include <fstream>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/repulsivegaussian.h"
#include "WaveFunctions/twoelectrons.h"
#include "WaveFunctions/manyelectrons.h"
#include "WaveFunctions/manyelectrons_coefficients.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorrepulsive.h"
#include "Hamiltonians/harmonicoscillatorelectrons.h"
#include "Hamiltonians/doubleharmonicoscillator.h"
#include "Hamiltonians/squarewell.h"
#include "Hamiltonians/finiteharmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "VariationMethods/steepestdescent.h"
#include "Math/random.h"
#include <mpi.h>
#include <cassert>

using namespace std;
using namespace UnitTest;


int main(int nargs, char* args[]) {

    int numprocs, my_rank;
    double totalE, totalKE, totalPE, totalVariance, totalAcceptanceRate, finalMeanDistance;
    double timeStart, timeEnd, totalTime;

    // Initialize MPI parallelization
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    timeStart = MPI_Wtime();

    int numberOfDimensions  = 2;
    int numberOfParticles   = 2;
    int numberOfSteps       = (int) 1e4;              // Monte Carlo cycles
    double omega            = 1.;                     // Oscillator frequency.
    double alpha            = 1.;//0.98456;//0.7;          // Variational parameter.         //3D: 0.983904
    double beta             = 0.40691;//2.82843;      // Variational parameter.         //3D: 0.376667
    double gamma            = 2.82843;
    double a                = 0.0043;                 // Hard core boson diameter.
    double stepLength       = 0.5;                    // Metropolis step length.
    double equilibration    = 0.1;                    // Amount of the total steps used for equilibration.
    double dt               = 0.01;                   // Time step for importance sampling.
    double aElectrons       = 1.; //1./3
    double C                = 1.;                     // Norm constant.
    double distToWall       = 3.;
    double V0               = 1.;

    vec L(3);
    L.fill(0.);
    L(0) = 5.;
    //L(1) = 5.;

    bool analyticalKinetic  = true;
    bool importanceSampling = true;
    bool repulsion          = false;                   // Switch for interacting system or not. (Coulomb for manybody qdot)
    bool quantumDots        = true;                   // Switch for quantum dot system.
    bool twobodyQD          = false;                  // Switch for twobody quantum dot system. (no Slater)
    bool Jastrow            = false;                   // Switch for Jastrow factor. (manybody qdot)
    bool optimizeParameters = false;                  // Switch for optimizing variational parameters.

    bool saveEnergies       = false;
    bool savePositions      = false;

    bool showProgress       = true;
    bool printToTerminal    = true;

    bool useCoeff 		    = false;				  // Coefficients c_ij = <ψ_i|φ_j> from the double well potential.

    bool doubleWell         = true;
    bool finiteWell         = false;
    bool squareWell         = false;

    bool runTests           = false;

    int numMyCycles = numberOfSteps/numprocs;

//    cout << "  -- Settings -- " << boolalpha << endl;
//    cout << " Analytical Kinetic : " << analyticalKinetic << endl;
//    cout << " Importance Sampling : " << importanceSampling << endl;
//    cout << " Repulsion : " << repulsion << endl;

    int numberOfDMCWalkers = 10; // Number of VMC configurations to save.
    std::vector<class Walker*> setOfWalkers(numberOfDMCWalkers);

    // Initiate System
    System* system = new System();

    // Set the set of walkers
    system->setWalkers(setOfWalkers);

    // RandomUniform creates a random initial state
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, my_rank));
    // Select which Hamiltonian and trial wave function to use (interacting or non-interacting)
    if (repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillatorRepulsive(system, omega, a, gamma, analyticalKinetic));
    system->setWaveFunction             (new RepulsiveGaussian(system, alpha, beta, a));
    }
    if (!repulsion && !quantumDots) {
    system->setHamiltonian              (new HarmonicOscillator(system, omega, analyticalKinetic));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    }
    if (quantumDots) {
        if (doubleWell) {
            system->setDoubleWellFlag   (doubleWell);
            system->setL                (L);
            system->setHamiltonian      (new DoubleHarmonicOscillator(system, L, omega, analyticalKinetic, repulsion));
        }
        else if (finiteWell) {
            system->setHamiltonian      (new FiniteHarmonicOscillator(system, distToWall, omega, analyticalKinetic, repulsion));
        }
        else if (squareWell) {
            system->setSquareWellFlag   (squareWell);
            system->setHamiltonian      (new SquareWell(system, V0, distToWall, omega, analyticalKinetic, repulsion));
        }
        else {
            system->setHamiltonian      (new HarmonicOscillatorElectrons(system, omega, analyticalKinetic, repulsion));
        }
        if (twobodyQD) {
            system->setWaveFunction     (new TwoElectrons(system, alpha, beta, omega, aElectrons, C, Jastrow));
        }
        else {
            if (useCoeff) {
                vec constants;
                system->retrieveConstantsFromFile("../diagonalization/PlotAndData/Constants.dat", constants);
                assert(omega == constants(0));
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
    // Optimize parameters
    if (optimizeParameters) {
        system->optimizeParameters          (system, alpha, beta);
    }
    system->setSaveEnergies             (saveEnergies);
    system->setSavePositions            (savePositions);
    // Start Monte Carlo simulation
    system->runMetropolisSteps          (numMyCycles, importanceSampling, showProgress, printToTerminal);

    // Compute MPI averages etc.
    system->MPI_CleanUp                 (totalE, totalKE, totalPE, totalVariance, totalAcceptanceRate, finalMeanDistance,
                                         timeStart, timeEnd, totalTime, numprocs, numberOfSteps);

    // Merge the files from the nodes into one data file
    system->mergeOutputFiles            (numprocs);


    //cout << setOfWalkers[0]->getParticles()[0]->getPosition()[0] << endl;
    //cout << setOfWalkers[1]->getParticles()[0]->getPosition()[0] << endl;

    if (runTests) { return RunAllTests(); }
    else { return 0; }
}

//12: 65.7
//20: 155.868

/*
 -- System info --
 Number of particles  : 50
 Number of dimensions : 3
 Number of Metropolis steps run : 10^6
 Number of equilibration steps  : 10^5

  -- Wave function parameters --
 Number of parameters : 2
 Parameter 1 : 0.5
 Parameter 2 : 2.82843

  -- Results --
 Energy : 127.306
 Variance : 0.191128
 Acceptance Rate : 0.891037

Computation Time : 13577

Optimizing alpha using steepest descent:

 -- System info --
 Number of particles  : 500
 Number of dimensions : 3
 Number of Metropolis steps run : 10^5
 Number of equilibration steps  : 10^4

  -- Wave function parameters --
 Number of parameters : 1
 Parameter 1 : 0.5

  -- Results --
 Energy : 750
 Variance : 7.64674e-06
 Acceptance Rate : 0.911399

Computation Time : 4580.39


*/
