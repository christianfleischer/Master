#ifndef PROJECT2_DMC_H
#define PROJECT2_DMC_H
#include <armadillo>

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
    int m_blockSize;


    class System* m_system;

    class Walker* m_trialWalker;
    std::vector<class Walker*> m_setOfWalkers;
};

#endif // PROJECT2_DMC_H
