#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "squarewell.h"

SquareWell::SquareWell(System* system, double V0, double distToWall, bool analyticalKinetic) : Hamiltonian(system, analyticalKinetic){
    m_distToWall = distToWall;
    m_V0 = V0;
}

std::vector<double> SquareWell::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double potentialEnergy = 0;

    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> r_i = particles[i]->getPosition();
        bool isOutside = false;
        for (int d = 0; d < numberOfDimensions; d++) {
            if (r_i[d] > m_distToWall) {
                isOutside = true;
            }
        }
        if (isOutside) potentialEnergy += m_V0;
    }

    double kineticEnergy = 0;

    if (m_analyticalKinetic == true){
        // Using analytical expression.
        double doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
        kineticEnergy = -0.5*doubleDerivative;
    }
    else{
        // Using numerical diff.
        kineticEnergy = computeKineticEnergy(particles);
    }

    std::vector<double> energies(3);
    energies[0] = kineticEnergy + potentialEnergy;
    energies[1] = kineticEnergy;
    energies[2] = potentialEnergy;

    return energies;
}
