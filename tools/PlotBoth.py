from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os


path = os.environ['DarkMatterDir']


d = 2


#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
    ax1.scatter(0,0,0, c='r')
elif  (d == 2):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0))
    ax1.scatter(0,0, c='r')




#Load data

RayData = glob.glob(path+'RT/*.txt')


def plot(x,y,z,PlotType):



#Plot it
    if (d == 3):

        if PlotType == 'plot':
            ax1.plot(x,y,z)  
        else:
            ax1.scatter(x,y,z, s=20)


        limit = max(max(abs(x)),max(abs(y)),max(abs(z)))
        ax1.set_xlim(-limit,+limit)
        ax1.set_ylim(-limit,+limit)
        ax1.set_zlim(-limit,+limit)


    if (d == 2):
        if PlotType == 'plot':
            ax1.plot(x,y)
        else:
            ax1.scatter(x,y,s=20)
     

        limit = max(max(abs(x)),max(abs(y)),max(abs(z)))
        ax1.set_xlim(-limit,+limit)
        ax1.set_ylim(-limit,+limit)


def PlotRay(file):
    
    data=np.loadtxt(file)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    plot(x,y,z,'plot')




    



for f in RayData:
    PlotRay(f)






#----------------------------

OrbitData = np.loadtxt(path + 'MPDFormatData.txt')


r = OrbitData[:,0]
theta = OrbitData[:,1]
phi = OrbitData[:,2]
a = OrbitData[0,3]

mm = np.sqrt(r**2 + a**2)
x = mm * np.sin(theta)*np.cos(phi)
y = mm * np.sin(theta)*np.sin(phi)
z = r  * np.cos(theta)


plot(x,y,z,'scatter')





plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

plt.show()
