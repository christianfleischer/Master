#include "doublewell.h"
#include "../system.h"

DoubleWell::DoubleWell(System *system, double omega)
    : WaveFunction(system, omega) {

}

void DoubleWell::harmonicOscillatorBasis(mat r, int n) {
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int N = m_system->getN();

    mat H = zeros(N-1, numberOfDimensions);
    for (int d = 0; d < numberOfDimensions; d++) {
        H.col(d) = computeHermitePolynomial(n, r.col(d));
    }
}

vec DoubleWell::potential (vec r, double L) {
    return m_omega*m_omega*(r%r - 2*abs(r)*L + L*L);
}
