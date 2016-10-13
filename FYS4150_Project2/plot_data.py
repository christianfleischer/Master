from pylab import *
import numpy as np
from math import factorial
import sys

# Position in first column(s), eigenfunction in the remainder:
table5nr = np.loadtxt('PlotAndData/omega5norepulsion.dat', skiprows=0)
table5constants = np.loadtxt('PlotAndData/omega5constants.dat', skiprows=0)
table5position = np.loadtxt('PlotAndData/omega5position.dat', skiprows=0)

omega, nDim, Lx, Ly, Lz, N, numEigFunctions = table5constants

numEigFunctions = int(numEigFunctions)
N = int(N)-1
nDim = int(nDim)

# Position vector:
r = np.zeros((nDim, N))
for d in range(nDim):
    r[d] = table5position[:,d]

# Numerically calculated eigenvector(2) (for different quantum numbers n):
psi = np.zeros((numEigFunctions,N))
for i in range(numEigFunctions):
    psi[i] = table5nr[:,i]


def H(r, n_r):
    factor = 2*r*np.sqrt(omega)
    HermitePolynomialPP = 0;                 
    HermitePolynomialP = 1;                  
    HermitePolynomial = HermitePolynomialP;  
    for n in range(1, n_r + 1):
        HermitePolynomial = factor*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    return HermitePolynomial

# Harmonic oscillator basis functions (for different quantum numbers n):
def phi(r, n):
    nFac = 1
    for i in n:
        nFac *= factorial(i)
    
    n2 = 2.**sum(n)
    
    pi4 = np.pi**(-0.25)
    const = omega**(0.25)*pi4/(np.sqrt(nFac*n2))
    rAbs2 = 0
    for i in r:
        rAbs2 += i*i

    wavefunc = np.exp(-0.5*omega*(rAbs2))
    

    phi = 1
    for d in range(nDim):
        phi *= H(r[d], n[d])

    return const*wavefunc*phi


# Number of quantum numbers to use:
nMax = 20

r[1] = np.zeros(N)
r[2] = np.zeros(N)

def C(r, n):
    return psi[0]*phi(r, n)

sup = 0
for n_x in range(nMax):
    sys.stdout.write('[%.1f %%]\r' % (n_x/nMax*100.))
    sys.stdout.flush()  
    for n_y in range(nMax):
        for n_z in range(nMax):
            n = [n_x, n_y, n_z]
            sup += C(r, n)*phi(r, n)

# When we set psi = psix*psiy*psiz, we must normalize manually:
print("<psi0|psi0>:   ", np.dot(psi[0],psi[0]))
print("<sup0|sup0>:   ", np.dot(sup,sup))

plot(r[0], sup**2/np.dot(sup,sup), r[0], psi[0]**2/np.dot(psi[0],psi[0])
, 'x')

title(r'''Probability density $|\psi(x)|^2$ with N=%d, $L_x = %.1f$,
        $L_, = %.1f$ and $x_{max/min}=\pm %d$''' 
        %(N+1, Lx, Ly, r[0,-1]+1))
ylabel(r'$|u(\rho)|^2$')
legend([r'Analytical solution', r'Numerical solution', r'Basis function solution'] 
        ,prop={'size':10})

show()
