#include "walker.h"
#include "dmc.h"
#include "../WaveFunctions/wavefunction.h"

Walker::Walker(int numberOfParticles, int numberOfDimensions)
{
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;

}

void Walker::setPosition(mat position) {
    m_r = position;
}

void Walker::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
    m_particles = initialState->getParticles();
}
