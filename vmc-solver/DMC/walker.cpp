#include "walker.h"
#include "dmc.h"

Walker::Walker(int numberOfParticles, int numberOfDimensions)
{
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;

}

void Walker::setPosition(mat position) {
    m_r = position;
}
