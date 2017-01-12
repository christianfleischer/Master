#ifndef DOUBLEHARMONICOSCILLATOR_H
#define DOUBLEHARMONICOSCILLATOR_H
#include "hamiltonian.h"
#include "../particle.h"
#include <vector>

class DoubleHarmonicOscillator : public Hamiltonian {
public:
    DoubleHarmonicOscillator(System* system, vec L, double omega, bool analyticalKinetic, bool repulsion);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);
    double evaluateSingleParticleWF(vec n, std::vector<double> r, int j);
    std::vector<double> computeSPWFDerivative(vec n, std::vector<double> r, int j);
    double computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j);
    double computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j);
    double computeHermitePolynomial(int nValue, double position);
    double computeHermitePolynomialDerivative(int nValue, double position);
    double computeHermitePolynomialDoubleDerivative(int nValue, double position);
    double computeHermitePolynomialAlphaDerivative(int nValue, double position);
    void setExpFactor(int randomParticle, std::vector<Particle *> particles);

private:
    int m_numberOfDimensions = 0;
    vec m_L;
    bool m_repulsion = false;

    class Hamiltonian* m_well1;
    class Hamiltonian* m_well2;

    std::vector<Particle *> m_particlesWell1;
    std::vector<Particle *> m_particlesWell2;
};

#endif // DOUBLEHARMONICOSCILLATOR_H
