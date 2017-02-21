#ifndef PROJECT2_HAMILTONIAN_H
#define PROJECT2_HAMILTONIAN_H
#include <vector>
#include <armadillo>

using namespace arma;

class Hamiltonian {
public:
    Hamiltonian(class System* system, bool analyticalKinetic);
    double computeKineticEnergy(std::vector<class Particle*> particles);
    virtual std::vector<double> computeLocalEnergy(std::vector<class Particle*> particles) { particles = particles; }
    virtual double evaluateSingleParticleWF(vec n, std::vector<double> r, int j) { return n[0]*r[0]*j; }
    virtual std::vector<double> computeSPWFDerivative(vec n, std::vector<double> r, int j) { j = j; n = n; return r; }
    virtual double computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j) { return n[0]*r[0]*j; }
    virtual double computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j) { return n[0]*r[0]*j; }
    virtual double computeHermitePolynomial(int n, double position ) {n = n; position = position;}
    double getAnalytic(){ return m_analyticalKinetic; }
    void setExpFactor(double expFactor) { m_expFactor = expFactor; }
    void setAlpha(double alpha) { m_alpha =  alpha; }
    virtual vec getWellDistance()   { return zeros(1); }

protected:
    class System* m_system = nullptr;
    bool m_analyticalKinetic = false;
    double m_expFactor = 0;
    double m_alpha = 0;
};

#endif // PROJECT2_HAMILTONIAN_H
