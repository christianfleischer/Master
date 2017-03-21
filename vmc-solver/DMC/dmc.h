#ifndef DMC_H
#define DMC_H
#include "walker.h"
#include "system.h"

class DMC
{
public:
    DMC(System* system, int N_c);
    void setEquilibrationSteps(int equilibration);
    void moveWalker(int iWalker);
    void runDMC();
private:
    int m_numberOfWalkers;
    int m_numberOfEquilibrationSteps = 0.1*m_numberOfWalkers;
    int m_numberOfParticles;
    int m_numberOfDimensions;

    System* m_system;

    Walker* m_trialWalker;
    Walker** m_setOfWalkers;
};

#endif // DMC_H
