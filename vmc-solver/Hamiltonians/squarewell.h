#ifndef SQUAREWELL_H
#define SQUAREWELL_H
#include "hamiltonian.h"
#include <vector>

class SquareWell : public Hamiltonian {
public:
    SquareWell(System* system, double V0, double distToWall, double omega, bool analyticalKinetic, bool repulsion);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);
    double evaluateSingleParticleWF(vec n, std::vector<double> r);
    std::vector<double> computeSPWFDerivative(vec n, std::vector<double> r);
    double computeSPWFDoubleDerivative(vec n, std::vector<double> r);
    double computeSPWFAlphaDerivative(vec n, std::vector<double> r);

private:
    double m_V0 = 0;
    double m_distToWall = 0;
    double m_omega = 0;
    bool m_repulsion = false;
};


#endif // SQUAREWELL_H
