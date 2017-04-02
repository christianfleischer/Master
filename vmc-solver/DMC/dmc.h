#ifndef PROJECT2_DMC_H
#define PROJECT2_DMC_H
#include <armadillo>

class DMC
{
public:
    DMC(class System* system, int N_c);
    void setEquilibrationSteps(double equilibration);
    void moveWalker(int iWalker);
    void runDMC();
    void copyWalker(class Walker* originalWalker, class Walker* newWalker);
    void copyParticles(std::vector<class Particle*> originalParticles, std::vector<class Particle*> newParticles);

    void setParameters(double alpha, double beta, double omega, double C, bool Jastrow);

private:
    int m_numberOfWalkers;
    int m_numberOfEquilibrationSteps = 0.1*m_numberOfWalkers;
    int m_numberOfParticles;
    int m_numberOfDimensions;
    int m_blockSize;

    double m_omega = 0;
    double m_alpha = 0;
    double m_beta = 0;
    double m_C = 0;

    bool m_Jastrow = false;

    class System* m_system;

    class Walker* m_trialWalker;
    std::vector<class Walker*> m_setOfWalkers;
};

#endif // PROJECT2_DMC_H
