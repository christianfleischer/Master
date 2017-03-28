#ifndef DMC_H
#define DMC_H
#include "walker.h"
#include "system.h"

class DMC
{
public:
    DMC(class System* system, int N_c);
    void setEquilibrationSteps(int equilibration);
    void moveWalker(int iWalker);
    void runDMC();
    void copyWalker(class Walker* originalWalker, Walker* newWalker);
private:
    int m_numberOfWalkers;
    int m_numberOfEquilibrationSteps = 0.1*m_numberOfWalkers;
    int m_numberOfParticles;
    int m_numberOfDimensions;


    class System* m_system;

    class Walker* m_trialWalker;
    std::vector<class Walker*> m_setOfWalkers;
};

#endif // DMC_H
