#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


void Potential (vec &V, vec r, double L, double omega)
{
    V = omega*omega*(r%r - 2*abs(r)*L + L*L);
    return;
}
