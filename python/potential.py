import numpy as np
import matplotlib.pylab as plt

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
