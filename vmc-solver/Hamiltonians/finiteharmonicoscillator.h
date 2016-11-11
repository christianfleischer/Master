#ifndef FINITEHARMONICOSCILLATOR_H
#define FINITEHARMONICOSCILLATOR_H
#include "hamiltonian.h"
#include <vector>


class FiniteHarmonicOscillator : public Hamiltonian {
public:
    FiniteHarmonicOscillator(System* system, double distToWall, double omega, bool analyticalKinetic, bool repulsion);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);
private:
    double m_omega = 0;
    double m_distToWall = 0;
    bool m_repulsion = false;
};

#endif // FINITEHARMONICOSCILLATOR_H
