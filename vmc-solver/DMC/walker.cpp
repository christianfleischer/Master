#include "walker.h"
#include "dmc.h"
#include "../WaveFunctions/wavefunction.h"

Walker::Walker(int numberOfParticles, int numberOfDimensions)
{
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;

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
