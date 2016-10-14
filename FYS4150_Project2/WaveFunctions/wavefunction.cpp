#include "wavefunction.h"

WaveFunction::WaveFunction(System* system, double omega) {
    m_system = system;
    m_omega = omega;
}

double WaveFunction::computeHermitePolynomial(int nValue, vec position) {
    // Computes Hermite polynomials.
    double omegaSqrt = sqrt(m_omega);
    double factor = 2*omegaSqrt;//*position;

    double HermitePolynomialPP = 0;                 // H_{n-2}
    double HermitePolynomialP = 1;                  // H_{n-1}
    double HermitePolynomial = HermitePolynomialP;  // H_n

    for (int n=1; n <= nValue; n++) {
        HermitePolynomial = factor*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    }

    return HermitePolynomial;
}
