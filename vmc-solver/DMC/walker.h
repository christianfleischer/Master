#ifndef WALKER_H
#define WALKER_H
#include <armadillo>
#include "system.h"
using namespace arma;


class Walker
{
public:
    Walker(int numberOfParticles, int numberOfDimensions);

    void killWalker() { m_isAlive = false; }
    void reviveWalker() { m_isAlive = true; }
    void setE(double energy) {m_E = energy; }

    double getE() { return m_E; }

    mat getPosition() { return m_r; };
    void setPosition(mat position);

    class System* getSystem()        { return m_system; }

private:
    bool m_isAlive = true;

    int m_numberOfParticles;
    int m_numberOfDimensions;

    double m_E;

    mat m_r;
    mat m_rAbs;

    class System* m_system = nullptr;

};

#endif // WALKER_H
