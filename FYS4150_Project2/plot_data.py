from pylab import *
import numpy as np
from math import factorial
import sys

# Position in first column(s), eigenfunction in the remainder:
table5nr = np.loadtxt('PlotAndData/omega5norepulsion.dat', skiprows=0)
table5constants = np.loadtxt('PlotAndData/omega5constants.dat', skiprows=0)
table5position = np.loadtxt('PlotAndData/omega5position.dat', skiprows=0)

omega, Lx, Ly, Lz, N, numEigFunctions = table5constants

numEigFunctions = int(numEigFunctions)
N = int(N)-1

# Position vector:
x,y,z = table5position[:,0], table5position[:,1], table5position[:,2]

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
def phi(x, n_x, y, n_y, z, n_z):
    nFac = factorial(n_x)*factorial(n_y)*factorial(n_z)
    n2 = 2.**(n_x + n_y + n_z)
    pi4 = np.pi**(-0.25)
    const = omega**(0.25)*pi4/(np.sqrt(nFac*2.**n_x))
    wavefunc = np.exp(-0.5*omega*(x*x + y*y + z*z))

    return const*H(x, n_x)*H(y, n_y)*H(z, n_z)*wavefunc


# Number of quantum numbers to use:
nMax = 20

# y = np.zeros(N)
# z = np.zeros(N)

def C(x, n_x, y, n_y, z, n_z):
    return psi[0]*phi(x, n_x, y, n_y, z, n_z)

sup = 0
for n_x in range(nMax):
    sys.stdout.write('[%.1f %%]\r' % (n_x/nMax*100.))
    sys.stdout.flush()  
    for n_y in range(nMax):
        for n_z in range(nMax):
            sup += C(x, n_x, y, n_y, z, n_z)*phi(x, n_x, y, n_y, z, n_z)

# When we set psi = psix*psiy*psiz, we must normalize manually:
print("<psi0|psi0>:   ", np.dot(psi[0],psi[0]))
print("<sup0|sup0>:   ", np.dot(sup,sup))

plot(x, sup**2/np.dot(sup,sup), x, psi[0]**2/np.dot(psi[0],psi[0])
, 'x')

title(r'''Probability density $|\psi(x)|^2$ with N=%d, $L_x = %.1f$,
        $L_, = %.1f$ and $x_{max/min}=\pm %d$''' 
        %(N+1, Lx, Ly, x[-1]+1))
ylabel(r'$|u(\rho)|^2$')
legend([r'Analytical solution', r'Numerical solution', r'Basis function solution'] 
        ,prop={'size':10})

show()
