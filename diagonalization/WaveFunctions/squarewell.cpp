#include "squarewell.h"
#include "../system.h"
#include "../Math/factorial.h"
#include <cmath>

SquareWell::SquareWell(System *system, double omega, double V0)
    : WaveFunction(system, omega) {

    m_V0 = V0;
}

vec SquareWell::harmonicOscillatorBasis(mat x, int n) {

    double nFac = factorial(n);

    double n2 = pow(2., n);
    double pi4 = pow(M_PI, -0.25);
    double omega4 = pow(m_omega, 0.25);
    double constant = omega4*pi4/sqrt(nFac*n2);

    vec xAbs2 = x%x;

    vec wavefunc = exp(-0.5*m_omega*xAbs2);

    vec phi = constant*wavefunc%computeHermitePolynomial(n, x);

    return phi;
}

vec SquareWell::potential (vec r, double L) {
    int N = m_system->getN();
    vec V(N+1);

    for (int i = 0; i < N+1; i++) {
        if (abs(r[i]) > L) { V[i] = m_V0; }
        else { V[i] = 0; }
    }

    return V;
}
