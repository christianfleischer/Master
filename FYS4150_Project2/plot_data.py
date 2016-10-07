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

# Numerically calculated eigenvector(2) (for different quantum numbers n):
psi = np.zeros((numEigFunctions,N))
for i in range(numEigFunctions-1):
    psi[i] = table5nr[:,i+1]

# Harmonic oscillator basis functions (for different quantum numbers n):
def phi(x, n_x):
    C = 12. # = m*w/hbar
    nFac = factorial(n_x)
    pi4 = np.pi**(-0.25)
    const = C**(0.25)*pi4/(np.sqrt(nFac*2**n_x))
    if n_x == 0:
        return const*np.exp(-C*0.5*x*x)
    if n_x == 1:
        return const*np.exp(-C*0.5*x*x) * C * 2*x
    if n_x == 2:
        return const*np.exp(-C*0.5*x*x) * C * (4*x*x - 2)
    if n_x == 3:
        return const*np.exp(-C*0.5*x*x) * C * (8*x*x*x - 12*x)
    if n_x == 4:
        return const*np.exp(-C*0.5*x*x) * C * (16*x*x*x*x - 48*x*x + 12)

# Coefficients for "analytical" super position psi_{i} = c_{ij}phi_{j}:
c = np.zeros(numEigFunctions - 1)

# Position of well minimum:
L = 2.5

# Super position of the basis functions (for different quantum numbers):
harOsc_plus0 = phi(x-L,0) + phi(x+L,0)
harOsc_minus0 = -phi(x-L,0) + phi(x+L,0)

harOsc_plus1 = phi(x-L,1) + phi(x+L,1)
harOsc_minus1 = -phi(x-L,1) + phi(x+L,1)

# Number of quantum numbers to use:
nMax = 4 
for n in range(0,nMax):
    c[n] = np.dot(psi[n], phi(x,n))/np.dot(phi(x,n), phi(x,n))
print(c)

print("<psi_0|psi_0>:   ", np.dot(psi[0],psi[0]))
print("<phi_0|phi_0>:   ", np.dot(phi(x,0),phi(x,0)))
print("<phi^+_0|phi^+_0>:   ", np.dot(harOsc_plus0,harOsc_plus0))


# Super position from the coefficient matrix c:
supPos = 0
for n in range(nMax):
    supPos += c[n]*phi(x,n)

plot(x, harOsc_plus0**2/np.dot(harOsc_plus0,harOsc_plus0), x, psi[0]**2)
# plot(x, (psi[0])**2, x, (supPos)**2)
# plot(x, phi(x,0),x, phi(x,1),x, phi(x,2),x, phi(x,3), x, supPos**2)
# hold('on')
# plot(x, harOsc_plus0**2, x, harOsc_minus0**2)

title(r'Probability density $|u(\rho)|^2$ for various $\omega_r$ with N=200 and $\rho_{max}$=|5|')
ylabel(r'$|u(\rho)|^2$')
legend([r'Analytical solution', r'Numerical solution', r'Basis function solution'] ,prop={'size':10})

show()
