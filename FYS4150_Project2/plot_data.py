from pylab import *
import numpy as np
from math import factorial
import sys

# Position in first column(s), eigenfunction in the remainder:

eigenvectors = np.loadtxt('PlotAndData/Eigenvectors.dat')

constants = np.loadtxt('PlotAndData/Constants.dat')
positionvectors = np.loadtxt('PlotAndData/Positionvectors.dat')

supCpp = np.loadtxt('PlotAndData/Superpositions.dat', skiprows=0)

omega, nDim, Lx, Ly, Lz, N, numEigFunctions, h = constants

numEigFunctions = int(numEigFunctions)
N = int(N)-1
nDim = int(nDim)

# Position vector:
r = np.zeros((nDim, N))
for d in range(nDim):
    r[d] = positionvectors[:,d]

# Numerically calculated eigenvector(2) (for different quantum numbers n):
psi = np.zeros((numEigFunctions,N))
for i in range(numEigFunctions):
    psi[i] = eigenvectors[:,i]


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
  

    pi4 = np.pi**(-0.25*nDim)
    const = omega**(0.25*nDim)*pi4/(np.sqrt(nFac*n2))
    rAbs2 = 0
    for i in r:
        rAbs2 += i*i

    wavefunc = np.exp(-0.5*omega*(rAbs2))
    

    phi = 1
    for d in range(nDim):
        phi *= H(r[d], n[d])

    return const*wavefunc*phi


# Number of quantum numbers to use:
nMax = 30

#r[1] = np.zeros(N)
#r[2] = np.zeros(N)

def C(r):
    nxprimeMax = 1
    C = np.zeros((nMax, nxprimeMax))

    for nxprime in range(nxprimeMax):
        for nx in range(nMax):
            innerprod = 0
            for i in range(N):
                innerprod += h*psi[nxprime,i]*phi(r,[nx])[i]
            C[nx,nxprime] = innerprod
    return C

def cOld(r,n):
    return psi[0]*phi(r, n)

#Cmat = C(r)

#a = []
#sup = 0
#for n_x in range(nMax):
#    sys.stdout.write('[%.1f %%]\r' % (float(n_x)/nMax*100.))
#    sys.stdout.flush()  
#    sup += Cmat[n_x,0]*phi(r, [n_x])

#    tmp1 = sup**2/np.dot(sup,sup) - psi[0]**2/np.dot(psi[0],psi[0])

#    a.append(np.sqrt(np.dot(tmp1, tmp1)))

# When we set psi = psix*psiy*psiz, we must normalize manually:
print("<psi0|psi0>:   ", np.dot(psi[0],psi[0]))
print("<supCpp0|supCpp0>:   ", np.dot(supCpp[:,0],supCpp[:,0]))

# plot(r[0], supCpp**2/np.dot(supCpp,supCpp), r[0], psi[0]**2/np.dot(psi[0],psi[0]), '+')


# plot(r[0], psi[0]**2/np.dot(psi[0],psi[0]), r[0], supCpp[:,0]**2/np.dot(supCpp[:,0],supCpp[:,0]), '+')
plot(r[0], psi[0]**2/np.dot(psi[0],psi[0]), r[0], supCpp[:,0]**2/np.dot(supCpp[:,0],supCpp[:,0]))

nVec = range(nMax)

# plot(nVec, a)


# hold('on')
# plot(r[0], psi[0]**2/np.dot(psi[0],psi[0]), r[0], sup**2/np.dot(sup,sup))
# plot(r[0], psi[0,:]**2, r[0], phi(r,[0]))

title(r'''Probability density $|\psi(x)|^2$ with N=%d, $L_x = %.1f$,
        $L_, = %.1f$ and $x_{max/min}=\pm %d$. ''' 
        %(N+1, Lx, Ly, r[0,-1]+1))
ylabel(r'$|u(\rho)|^2$')
legend([r'Analytical solution', r'Numerical solution', r'Basis function solution'] 
        ,prop={'size':10})

show()
