#ifndef FINITEWELL_H
#define FINITEWELL_H

#include "wavefunction.h"

class FiniteWell : public WaveFunction {
public:
    FiniteWell(class System* system, double omega);
    vec harmonicOscillatorBasis(mat x, int n);
    vec potential (vec r, double L);
};

#endif // FINITEWELL_H
