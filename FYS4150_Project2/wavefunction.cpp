#include "wavefunction.h"

waveFunction::waveFunction() {

}


//void waveFunction::harmonicOscillatorBasis(){
//    mat H = zeros(N-1, m_numberOfDimensions);
//    for (int d = 0; d < m_numberOfDimensions; d++) {
//        H.col(d) = computeHermitePolynomial(nValue, r.col(d));
//    }
//}
//
//double waveFunction::computeHermitePolynomial(int nValue, double position) {
//    // Computes Hermite polynomials.
//    double omegaSqrt = sqrt(m_omega);
//    double factor = 2*omegaSqrt*position;
//
//    double HermitePolynomialPP = 0;                 // H_{n-2}
//    double HermitePolynomialP = 1;                  // H_{n-1}
//    double HermitePolynomial = HermitePolynomialP;  // H_n
//
//    for (int n=1; n <= nValue; n++) {
//        HermitePolynomial = factor*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
//        HermitePolynomialPP = HermitePolynomialP;
//        HermitePolynomialP = HermitePolynomial;
//    }
//
//    return HermitePolynomial;
//}
