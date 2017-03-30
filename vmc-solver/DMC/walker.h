#ifndef PROJECT2_WALKER_H
#define PROJECT2_WALKER_H
#include <armadillo>
#include "system.h"
#include "dmc.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../WaveFunctions/wavefunction.h"

using namespace arma;


class Walker
{
public:
    Walker(int numberOfParticles, int numberOfDimensions);

    void                            killWalker() { m_isAlive = false; }
    void                            reviveWalker() { m_isAlive = true; }
    void                            setE(double energy) {m_E = energy; }
    void                            setPosition(mat position);
    void                            setInitialState(InitialState* initialState);

    int                             getNumberOfParticles() { return m_numberOfParticles; }
    int                             getNumberOfDimensions() { return m_numberOfDimensions; }
    double                          getE() { return m_E; }
    mat                             getPosition() { return m_r; };
    class InitialState*             getInitialState()   { return m_initialState; }
    class System*                   getSystem()        { return m_system; }
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }

private:
    bool                            m_isAlive = true;
    bool                            m_saveEnergies = false;
    bool                            m_savePositions = false;

    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_currentParticle = 0;

    double                          m_E = 0;
    double                          m_stepLength = 0.1;
    double                          m_dt = 0.01;
    double                          m_computationTime = 0;

    mat                             m_r;
    mat                             m_rAbs;

    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();

    class WaveFunction*             m_waveFunction = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    class System*                   m_system = nullptr;

    FILE*                           m_outfileE;
    FILE*                           m_outfileP;
};

#endif // PROJECT2_WALKER_H
