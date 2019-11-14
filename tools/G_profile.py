from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os

#Universal constants
Msolar = 1.989e30
G = 6.67e-11
c=3e8
MBH = 4.31e6
pc = 3.0860e16


#Conversion parameter
convert = c**2/(G*MBH*Msolar)

# DM consants
r0_SI = 10 #kpc
r0 = r0_SI * 1000 *pc#* convert
rho0 = 4e7
k = rho0 * r0_SI**3

print ('k = ', k)



#Define G

r = np.linspace(1000*pc,8000*pc,100)
RR = r/r0
kappa = -2*k*np.pi/r
A = 1 + RR



G = (1+RR**2)**(kappa*(1-RR))*A**(2*kappa*A)*np.exp(-2*kappa*A*np.arctan(RR)) - 2/r

Y = r * (1-G) / 2



r1 = 0.1
r = np.linspace(0.01,50,100)
den = 10*r1/((r+r1)*(r**2+r1**2))


#aa = 0.2**2
#Y = r**2 * G + aa

plt.figure()
plt.plot(r,den)
#plt.ylim(-1,0.5)
plt.show()








