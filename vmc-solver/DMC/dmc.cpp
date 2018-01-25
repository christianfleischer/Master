#include "dmc.h"
#include "walker.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/manyelectronsDMC.h"
#include "system.h"
#include "Hamiltonians/hamiltonian.h"
#include "particle.h"
#include <armadillo>

DMC::DMC(System* system, int N_c)
{
    m_numberOfWalkers = N_c;
    m_system = system;
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles = m_system->getNumberOfParticles();

    m_trialWalker = new Walker(m_numberOfParticles, m_numberOfDimensions);

    m_setOfWalkers = m_system->getWalkers();

}

void DMC::runDMC() {

    //equlibration steps:
    for (int cycle = 1; cycle <= m_numberOfEquilibrationSteps; cycle++) {
        for (int iWalker = 0; iWalker < m_numberOfWalkers; iWalker++) {
            moveWalker(iWalker);
        }
    }

    //actual dmc steps:
    for (int cycle = 1; cycle <= m_numberOfWalkers; cycle++) {
        for (int iWalker = 0; iWalker < m_numberOfWalkers; iWalker++) {
            moveWalker(iWalker);
        }

    }


}

void DMC::setEquilibrationSteps(double equilibration) {
    m_numberOfEquilibrationSteps = equilibration*m_numberOfWalkers;
}

void DMC::copyWalker(Walker* originalWalker, Walker* newWalker) {
    copyParticles(originalWalker->getParticles(), newWalker->getParticles());

    //Copy wavefunction:
    newWalker->setWaveFunction(new ManyElectronsDMC(m_system, m_alpha, m_beta, m_omega, m_C, m_Jastrow));

    newWalker->setLocalE(originalWalker->getLocalE());

    //These are the values that are needed in the wavefunction:
    mat SPWFMatOld = originalWalker->getWaveFunction()->getSPWFMat();
    field<vec> SPWFDMatOld = originalWalker->getWaveFunction()->getSPWFDMat();
    mat SPWFDDMatOld = originalWalker->getWaveFunction()->getSPWFDDMat();

    newWalker->getWaveFunction()->setSPWFMat(SPWFMatOld);
    newWalker->getWaveFunction()->setSPWFDMat(SPWFDMatOld);
    newWalker->getWaveFunction()->setSPWFDDMat(SPWFDDMatOld);
}

double getBranchingFunction(double E, double EOld, double E_T)
{
    return exp(-(0.5*(E + EOld) - E_T));
}

void DMC::moveWalker(int currentWalker) {
    copyWalker(m_setOfWalkers[currentWalker], m_trialWalker);

    for (int block = 0; block < m_blockSize; block++) {
        double localEOld = m_setOfWalkers[currentWalker]->getLocalE();
        std::vector<double> energies = m_system->getHamiltonian()->computeLocalEnergy(m_setOfWalkers[currentWalker]->getParticles());
        double localEnergy = energies[0];
        m_setOfWalkers[currentWalker]->setLocalE(localEnergy);

        double E_T = 0;
        double GB = getBranchingFunction(localEnergy, localEOld, E_T);

    }
}


void DMC::copyParticles(std::vector<Particle *> originalParticles, std::vector<Particle *> newParticles) {
    for (int i = 0; i < m_numberOfParticles; i++) {
        newParticles[i]->setPosition(originalParticles[i]->getPosition());
    }
}


void DMC::setParameters(double alpha, double beta, double omega, double C, bool Jastrow) {
    m_alpha = alpha;
    m_beta = beta;
    m_omega = omega;
    m_C = C;
    m_Jastrow = Jastrow;
}

