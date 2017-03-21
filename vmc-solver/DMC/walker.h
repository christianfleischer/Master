#ifndef WALKER_H
#define WALKER_H


class Walker
{
public:
    Walker(int numberOfParticles, int numberOfDimensions);

    void killWalker() { m_isAlive = false; }
    void reviveWalker() { m_isAlive = true; }
    void setE(double energy) {m_E = energy; }

    double getE() { return m_E; }


private:
    bool m_isAlive = true;

    int m_numberOfParticles;
    int m_numberOfDimensions;

    double m_E;
};

#endif // WALKER_H
