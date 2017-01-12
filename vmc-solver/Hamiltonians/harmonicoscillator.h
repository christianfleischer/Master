#ifndef PROJECT2_HARMONICOSCILLATOR_H
#define PROJECT2_HARMONICOSCILLATOR_H
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, bool analyticalKinetic);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);

private:

};

#endif // PROJECT2_HARMONICOSCILLATOR_H
