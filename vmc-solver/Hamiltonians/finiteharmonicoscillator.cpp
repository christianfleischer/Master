#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "finiteharmonicoscillator.h"

FiniteHarmonicOscillator::FiniteHarmonicOscillator(System *system, double distToWall, double omega, bool analyticalKinetic, bool repulsion) : Hamiltonian(system, analyticalKinetic){
    assert(omega > 0);
    m_omega = omega;
    m_repulsion = repulsion;
    m_distToWall = distToWall;
}

std::vector<double> FiniteHarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
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
    double repulsiveTerm = 0;

    for (int i=0; i < numberOfParticles; i++){
        double rSquared = 0;
        std::vector<double> r_i = particles[i]->getPosition();
        for (int k=0; k < numberOfDimensions; k++){
            rSquared += r_i[k]*r_i[k];
        }
        if (rSquared > m_distToWall*m_distToWall) { rSquared = m_distToWall*m_distToWall; }
        potentialEnergy += rSquared;

        for (int j=i+1; j < numberOfParticles; j++){
            double r_ijSquared = 0;
            std::vector<double> r_j = particles[j]->getPosition();
            for (int k=0; k < numberOfDimensions; k++){
                    r_ijSquared += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
            }

            double r_ij = sqrt(r_ijSquared);
            repulsiveTerm += 1./r_ij;
        }
    }
    potentialEnergy *= 0.5*m_omega*m_omega;
    if (m_repulsion) { potentialEnergy += repulsiveTerm; }

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

