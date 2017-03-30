#include "dmc.h"
#include "walker.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/manyelectronsDMC.h"
#include "system.h"
#include <armadillo>

DMC::DMC(System* system, int N_c)
{
    m_numberOfWalkers = N_c;
    m_system = system;

    m_trialWalker = new Walker(m_numberOfParticles, m_numberOfDimensions);
    m_setOfWalkers = m_system->getWalkers();

    for (int i = 0; i < m_numberOfWalkers; i++) {
        m_setOfWalkers[i] = new Walker(m_numberOfParticles, m_numberOfDimensions);
    }
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

void DMC::setEquilibrationSteps(int equilibration) {
    m_numberOfEquilibrationSteps = equilibration;
}

void DMC::copyWalker(Walker* originalWalker, Walker* newWalker) {
    newWalker->setPosition(originalWalker->getPosition());
    newWalker->setE(originalWalker->getE());

    //These are the values that are needed in the wavefunction:
    mat SPWFMatOld = originalWalker->getWaveFunction()->getSPWFMat();
    field<vec> SPWFDMatOld = originalWalker->getWaveFunction()->getSPWFDMat();
    mat SPWFDDMatOld = originalWalker->getWaveFunction()->getSPWFDDMat();

    newWalker->getWaveFunction()->setSPWFMat(SPWFMatOld);
    newWalker->getWaveFunction()->setSPWFDMat(SPWFDMatOld);
    newWalker->getWaveFunction()->setSPWFDDMat(SPWFDDMatOld);
}

void DMC::moveWalker(int currentWalker) {
    copyWalker(m_setOfWalkers[currentWalker], m_trialWalker);
    for (int block = 0; block < m_blockSize; block++) {
        double localE = m_setOfWalkers[currentWalker]->getE();

        for (int currentParticle = 0; currentParticle < m_numberOfParticles; currentParticle++){
            bool acceptedStep;
            // Ask Metropolis step functions whether the step was accepted or not
            acceptedStep = m_system->metropolisStepImpSamplingDMC(currentParticle, m_trialWalker);

            //This is the sampling part of the iteration:

            //if (i >= m_trialWalker->getSystem()->getNumberOfMetropolisSteps()){
            //    // Sample energy etc.
            //    m_trialWalker->getSystem()->m_sampler->sample(acceptedStep);
            //}
        }
    }
}
