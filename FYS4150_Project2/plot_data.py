from pylab import *
import numpy as np
from math import factorial

# Position in first column(s), eigenfunction in the remainder:
table5nr = np.loadtxt('PlotAndData/omega5norepulsion.dat', skiprows=0)

# Number of x values:
N = len(table5nr[:,0])

# Number of eigenfunctions fetched from the solver:
numEigFunctions = len(table5nr[0,:])

# Position vector:
x = table5nr[:,0]
print(x[0])

# Numerically calculated eigenvector(2) (for different quantum numbers n):
psi = np.zeros((numEigFunctions,N))
for i in range(numEigFunctions-1):
    psi[i] = table5nr[:,i+1]

# Harmonic oscillator basis functions (for different quantum numbers n):
def phi(x, n_x):
    nFac = factorial(n_x)
    pi4 = np.pi**(-0.25)
    const = pi4/(np.sqrt(nFac*2**n_x))
    if n_x == 0:
        return const*np.exp(-0.5*x*x)
    if n_x == 1:
        return const*np.exp(-0.5*x*x) * x

# Coefficients for "analytical" super position psi_{i} = c_{ij}phi_{j}:
c = np.zeros(numEigFunctions - 1)

# Position of well minimum:
L = 5.

# Super position of the basis functions (for different quantum numbers):
harOsc_plus0 = phi(x-L,0) + phi(x+L,0)
harOsc_minus0 = -phi(x-L,0) + phi(x+L,0)

harOsc_plus1 = phi(x-L,1) + phi(x+L,1)
harOsc_minus1 = -phi(x-L,1) + phi(x+L,1)


# Number of quantum numbers to use:
nMax = 2 

for i in range(0,nMax):
    c[i] = np.dot(psi[i], phi(x,i))
print(c)

# Super position from the coefficient matric c:
supPos = 0
for n in range(2):
    supPos += c[n]*phi(x,n)

# plot(x, harOsc_plus0, x, harOsc_minus0)
# plot(x, (psi[0])**2, x, (supPos)**2)
plot(x, psi[2]**2, x, psi[3]**2)
hold('on')
# plot(x, harOsc_plus0**2, x, harOsc_minus0**2)

title(r'Probability density $|u(\rho)|^2$ for various $\omega_r$ with N=200 and $\rho_{max}$=|5|')
ylabel(r'$|u(\rho)|^2$')
legend([r'Analytical solution', r'Numerical solution', r'Basis function solution'] ,prop={'size':10})

show()
