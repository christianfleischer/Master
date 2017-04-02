#include "walker.h"
#include "dmc.h"
#include "../WaveFunctions/wavefunction.h"
#include "particle.h"

Walker::Walker(int numberOfParticles, int numberOfDimensions)
{
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    for (int i = 0; i < m_numberOfParticles; i++) {
        m_particles.push_back(new Particle());
        m_particles[i]->setNumberOfDimensions(m_numberOfDimensions);
    }
}

void Walker::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
    m_particles = initialState->getParticles();
}

void Walker::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void Walker::setParticles(std::vector<Particle *> particles) {
    m_particles = particles;
}

void Walker::setSystem(System *system) {
    m_system = system;
}

