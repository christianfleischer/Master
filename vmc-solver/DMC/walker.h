#ifndef PROJECT2_WALKER_H
#define PROJECT2_WALKER_H
#include <armadillo>
#include "../InitialStates/initialstate.h"
#include "system.h"

using namespace arma;


class Walker
{
public:
    Walker(int numberOfParticles, int numberOfDimensions);

    void                            killWalker() { m_isAlive = false; }
    void                            reviveWalker() { m_isAlive = true; }
    void                            setLocalE(double energy) { m_E = energy; }
    void                            setInitialState(InitialState* initialState);
    void                            setWaveFunction(class WaveFunction* waveFunction);
    void                            setParticles(std::vector<class Particle*> particles);
    void                            setSystem(class System* system);

    class Sampler*                  m_sampler = nullptr;

    int                             getNumberOfParticles() { return m_numberOfParticles; }
    int                             getNumberOfDimensions() { return m_numberOfDimensions; }
    double                          getLocalE() { return m_E; }
    class InitialState*             getInitialState()   { return m_initialState; }
    class System*                   getSystem()        { return m_system; }
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class WaveFunction*             m_waveFunction = nullptr;

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

    std::vector<class Particle*>    m_particles;// = std::vector<class Particle*>();

    class InitialState*             m_initialState = nullptr;
    class System*                   m_system = nullptr;

    FILE*                           m_outfileE;
    FILE*                           m_outfileP;
};

#endif // PROJECT2_WALKER_H
