import numpy as np
import matplotlib.pylab as plt
import sys
'''
N = 200 

h = 10./N

L = [2.5,0,0]
o = 1.

rho = -5 + np.linspace(0, N, N+1)*h;
V = 0.5*o*o*(np.square(rho)-2*abs(rho)*L[0]+L[0]*L[0]);

for i in range(int(N/2)):
    V[i] += 2;

plt.plot(rho,V)
plt.show()
'''
positionvectors = np.loadtxt('../diagonalization/PlotAndData/Positionvectors.dat')
potential = np.loadtxt('../diagonalization/PlotAndData/Potential.dat')
constants = np.loadtxt('../diagonalization/PlotAndData/Constants.dat')

omega, nDim, Lx, Ly, Lz, N, numEigFunctions, h = constants
N = int(N)-1

# Position vector:
r = np.zeros((nDim, N))
for d in range(int(nDim)):
    r[d] = positionvectors[:,d]

plt.plot(r[0], potential[1:-1])
plt.show()

for i in range(N):
    print potential[i,0]




