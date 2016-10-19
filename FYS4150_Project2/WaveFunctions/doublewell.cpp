#include "doublewell.h"
#include "../system.h"
#include "../Math/factorial.h"
#include <cmath>

DoubleWell::DoubleWell(System *system, double omega)
    : WaveFunction(system, omega) {

}

vec DoubleWell::harmonicOscillatorBasis(mat r, vec n) {
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int N = m_system->getN();

    double nFac = 1.;
    int nSum = 0;

    for (int d = 0; d < numberOfDimensions; d++) {
        nFac *= factorial(n(d));
        nSum += n(d);
    }

    double n2 = pow(2., nSum);
    double pi4 = pow(M_PI, -0.75);
    double omega4 = pow(m_omega, 0.75);
    double constant = omega4*pi4/sqrt(nFac*n2);

    vec rAbs2 = zeros(N-1);

    for (int d = 0; d < numberOfDimensions; d++) {
        rAbs2 += r.col(d)%r.col(d);
    }

    vec wavefunc = exp(-0.5*m_omega*rAbs2);

    vec phi = ones(N-1);
    for (int d = 0; d < numberOfDimensions; d++) {
        phi %= computeHermitePolynomial(n(d), r.col(d));
    }

    return constant*wavefunc%phi;
}

vec DoubleWell::potential (vec r, double L) {
    return m_omega*m_omega*(r%r - 2*abs(r)*L + L*L);
}
