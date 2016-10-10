from pylab import *
import numpy as np
from math import factorial

# Position in first column(s), eigenfunction in the remainder:
table5nr = np.loadtxt('PlotAndData/omega5norepulsion.dat', skiprows=0)
table5constants = np.loadtxt('PlotAndData/omega5constants.dat', skiprows=0)
table5position = np.loadtxt('PlotAndData/omega5position.dat', skiprows=0)

omega, Lx, Ly, N, numEigFunctions = table5constants

numEigFunctions = int(numEigFunctions)
N = int(N)-1

# Position vector:
x,y = table5position[:,0], table5position[:,1]

# Numerically calculated eigenvector(2) (for different quantum numbers n):
psi = np.zeros((numEigFunctions,N))
for i in range(numEigFunctions):
    psi[i] = table5nr[:,i]

# Harmonic oscillator basis functions (for different quantum numbers n):
omega = 5.
def phi(x, n_x):
    nFac = factorial(n_x)
    pi4 = np.pi**(-0.25)
    const = omega**(0.25)*pi4/(np.sqrt(nFac*2.**n_x))
    factor = 2*x*np.sqrt(omega)

    HermitePolynomialPP = 0;                 
    HermitePolynomialP = 1;                  
    HermitePolynomial = HermitePolynomialP;  
    for n in range(1, n_x + 1):
        HermitePolynomial = factor*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    return const*HermitePolynomial*np.exp(-0.5*omega*x*x)

# Coefficients for "analytical" super position psi_{i} = c_{ij}phi_{j}:
c = np.zeros(numEigFunctions - 1)


# Number of quantum numbers to use:
nMax = numEigFunctions-1

def C(x,n_x):
    return psi[0]*phi(x, n_x)

print("<psi_0|psi_0>:   ", np.dot(psi[0],psi[0]))
print("<phi_0|phi_0>:   ", np.dot(phi(x,0),phi(x,0)))

n_x = 0
# plot(x, harOsc_plus(x,n_x)**2/np.dot(harOsc_plus(x,n_x),harOsc_plus(x,n_x)),
#         'x', x, psi[2*n_x]**2)

sup = 0
for i in range(100):
    sup += C(x,i)*phi(x,i)

print("<sup0|sup0>:   ", np.dot(sup,sup))

plot(x, sup**2/np.dot(sup,sup), x, psi[0]**2)

title(r'''Probability density $|\psi(x)|^2$ with N=%d, $L_x = %.1f$,
        $L_, = %.1f$, $n_x = %d$ and $x_{max/min}=\pm %d$''' 
        %(N+1, Lx, Ly, n_x, x[-1]+1))
ylabel(r'$|u(\rho)|^2$')
legend([r'Analytical solution', r'Numerical solution', r'Basis function solution'] 
        ,prop={'size':10})

show()
