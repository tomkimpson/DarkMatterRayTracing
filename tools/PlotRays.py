from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os


path = os.environ['RayTracingDir']


d = 2


#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0))



#Load data
data = np.loadtxt(path+'test.txt')

x = data[:,0]
y = data[:,1]
z = data[:,2]



#Plot it


if (d == 3):
    ax1.plot(x,y,z)  
    limit = max(max(x),max(y),max(z))
    ax1.set_xlim(-limit,+limit)
    ax1.set_ylim(-limit,+limit)
    ax1.set_zlim(-limit,+limit)


    #Plot BH
    ax1.scatter(0,0,0, c='r')

    #Plot PSR
    ax1.scatter(x[0],y[0],z[0], c='b')






if (d == 2):
    ax1.plot(x,y)

    limit = max(max(x),max(y),max(z))
    ax1.set_xlim(-limit,+limit)
    ax1.set_ylim(-limit,+limit)


    #Plot BH
    ax1.scatter(0,0, c='r')

    #Plot PSR
    ax1.scatter(x[0],y[0], c='b')




plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

plt.show()
