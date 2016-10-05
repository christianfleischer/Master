from pylab import *
import numpy as np

table5nr = np.loadtxt('PlotAndData/omega5norepulsion.dat', skiprows=0)
x = table5nr[:,0]
N = len(table5nr[:,0])
numEn = len(table5nr[0,:])
print(numEn)

u5nr = np.zeros((numEn,N))

for i in range(numEn-1):
    u5nr[i] = table5nr[:,i+1]

omega_r = 5

def phi(x, n_x):
    if n_x == 0:
        return np.exp(-0.5*omega_r*x*x)
    if n_x == 1:
        return 2*np.sqrt(omega_r)*x*np.exp(-0.5*omega_r*x*x)


c = np.zeros(numEn - 1)

harOsc_plus1 = phi(x-2.5,1) + phi(x+2.5,1)
harOsc_minus1 = -phi(x-2.5,1) + phi(x+2.5,1)

for n in range(1,3):
    c[n-1] = np.dot(table5nr[:,n], phi(x,n-1))
    
print(c)

supPos = 0
for n in range(2):
    supPos += c[n]*phi(x,n)

plot(x, u5nr[2], x, harOsc_plus1)
title(r'Probability density $|u(\rho)|^2$ for various $\omega_r$ with N=200 and $\rho_{max}$=|5|')
ylabel(r'$|u(\rho)|^2$')
legend([r'$\omega_r$=0.01, No repulsion'] ,prop={'size':10})
show()
