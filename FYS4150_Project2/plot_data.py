from pylab import *
import numpy as np

table0_01nr = np.loadtxt('PlotAndData/omega0_01norepulsion.dat', skiprows=0)

x = table0_01nr[:,0]
u0_01nr = table0_01nr[:,1]

plot(x, u0_01nr**2)
title(r'Probability density $|u(\rho)|^2$ for various $\omega_r$ with N=200 and $\rho_{max}$=|5|')
ylabel(r'$|u(\rho)|^2$')
legend([r'$\omega_r$=0.01, No repulsion'] ,prop={'size':10})

show()
