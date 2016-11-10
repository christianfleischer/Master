#ifndef SQUAREWELL_H
#define SQUAREWELL_H
#include "hamiltonian.h"
#include <vector>

class SquareWell : public Hamiltonian {
public:
    SquareWell(System* system, double V0, double distToWall, bool analyticalKinetic);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);
private:
    double m_V0 = 0;
    double m_distToWall = 0;
};


#endif // SQUAREWELL_H
