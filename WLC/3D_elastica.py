#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 13:41:44 2021

@author: michaelselby
"""

import scipy
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tikzplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from matplotlib import cm




# parameters
# KT is kBT in nmpN
# a1 is the bending modulus of the floppy segment with units 
# L is the length of the floppy segment in nm
# a is the bending modulus of the rest of the curve
# f is the applied force in pN


KT = 4.28
A = 45*KT
C = 100*KT
f = 2
bp = 100
L = bp*0.34  #length in nm
w = 0.5
npoints=1000
wpoints = 2
Lk = -10
p = 3.4
Lk0 = bp/p
lam = 2
u0 = 0.9
phi0 = 0
psi0 = 0
 
# define the 4 Jacobi elliptic functions: sn, cn, dn and am
def sn(u,m):
    return scipy.special.ellipj(u,m)[0]
def cn(u,m):
    return scipy.special.ellipj(u,m)[1]
def dn(u,m):
    return scipy.special.ellipj(u,m)[2]
def am(u,m):
    return scipy.special.ellipj(u,m)[3]
# define the incomplete elliptic integral of the second kind
def E(theta,m):
    return scipy.special.ellipeinc(theta,m)
def F(theta,m):
    return scipy.special.ellipkinc(theta,m)

svals = np.linspace(-L/2,L/2,npoints)
svalspos = np.linspace(0,L/2,npoints)
wvals = np.linspace(-w/2,w/2,wpoints)

ds = L/(npoints-1)






'''
Solving the Restraints

The type of constraints depend on what type of solution we are searching for

For closed loops, f is a lagrenge multiplier which must be determined, whereas for MT loops, it is an input

'''



##                      Figure of 8                          ##
###############################################################

def mfig8(a,b,c):
    return (b-c)/(a-c)

def fig8theta(a,b,c,f,A,s):
    lamd = (A/f)**0.5
    return np.arccos((b-c)*(sn(s*(a-c)**0.5/(2*lamd)**0.5+L/(4*lamd),mfig8(a,b,c))+c))





##                      Loopus Baybee                        ##
###############################################################


def psi(A,C,L,Lk):
    return A*C*(Lk+1)/((A+C)*L)

def function(b,A,L,f,psi):
    u = L*(1+psi**2/(2*(1-b)))**0.5/(2*lam)
    return (b+1)*(sn(u,m(b,psi)))**2-2

def m(b,psi):
    if (psi<0.001):
        return (b+1)/2
    else:
        return (1-b**2)/(psi**2+2*(1-b))

def theta_loop(s,b,m,A,L,f,psi):
    u = s*((b+1)/2)**0.5/(lam*m**0.5)
    if s>=0:
        return np.arccos((b+1)*(sn(u,m))**2-1)
    else:
        return -np.arccos((b+1)*(sn(u,m))**2-1)

def phi_loop(s,b,m,A,L,f,psi):
    u = (s*((b+1)/2)**0.5)/(lam*m**0.5)
    return phi0-psi*(2*m**0.5)**0.5/((A*f)**0.5*(b+1)**1.5)*(s-cn(u,m)*dn(u,m)/sn(u,m)-E(am(u,m),m))


def theta_homoclinic(s,b,m,A,L,f,psi):
    u = s/lam
    return np.arccos(2*(sn(u,m))**2-1)


psi=psi(A,C,L,Lk)



b_exact = fsolve(lambda b: function(b,A,L,f,psi),0.1)[0]
# b_exact = 1
m_exact = m(b_exact,psi)
# m_exact = 1


'''
Now we construct the ribbon 
'''

# define the rotation matrix R1
def rot1(theta):
    return np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])

# defines the rotation matrix R3
def rot3(theta):
    return np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])


# defines the rotation R3(psi)*R1(theta)*R3(phi)
def rot(psi,theta,phi):
    return np.matmul(np.matmul(rot3(psi),rot1(theta)),rot3(phi))


# initial tangent vectors
e30 = np.array([0,0,1])
e20 = np.array([0,1,0])
e10 = np.array([1,0,0])



psivals = []
thetavals = []
phivals = []
for s in svals:
    thetavals.append(theta_loop(s,b_exact, m_exact, A, L, f, psi))
    # thetavals.append(theta_homoclinic(s,b_exact, m_exact, A, L, f, psi))
    phivals.append(phi_loop(s,b_exact, m_exact, A, L, f, psi))
    psivals.append(psi0+psi*(1/A+1/C)*s-phi_loop(s,b_exact, m_exact, A, L, f, psi))

    



tangentvals = []
normalvals = []
binormalvals = []
for i in range(0,len(svals)):
    tangentvals.append(np.dot(rot(psivals[i],thetavals[i],phivals[i]),e30))
    normalvals.append(np.dot(rot(psivals[i],thetavals[i],phivals[i]),e20))
    binormalvals.append(np.dot(rot(psivals[i],thetavals[i],phivals[i]),e10))



r = np.array([0,0,-L/2])
rvals = []
for i in range(0,len(svals)):
    if (i==0):
        rvals.append(r)
    else:
        r = r + tangentvals[i]*ds   # position is integral of the tangent
        rvals.append(r)

zcom = 0
for i in range(0,npoints):
    zcom = zcom + rvals[i][2]
zcom = zcom/npoints
               
             
print("k = " +str(m_exact**0.5))
print("a = " +str((b_exact+1)/m_exact-1))
print("b = " +str(b_exact))



'''
Now for the plot

'''


fig = plt.figure()
ax = fig.gca(projection='3d')

xvals = []
yvals = []
zvals = []
for j in range(0,len(rvals)):
        xvals.append(rvals[j][0])  # recenter loop
        yvals.append(rvals[j][1])      
        zvals.append(rvals[j][2]-zcom)
        
xsep  = max(xvals)-min(xvals)
ysep  = max(yvals)-min(yvals)
zsep  = max(zvals)-min(zvals)

sep = max([xsep,ysep,zsep])

xvalslong = []
yvalslong = []
zvalslong = []

for j in range(0,len(rvals)):
    for i in range(0,wpoints):
        xvalslong.append(rvals[j][0]+wvals[i]*binormalvals[i][0])
        yvalslong.append(rvals[j][1]+wvals[i]*binormalvals[i][1])
        zvalslong.append(rvals[j][2]+wvals[i]*binormalvals[i][2]-zcom)   
        
ax.scatter(zvals,xvals,yvals, c='red',s=0.1)
ax.plot3D(zvalslong,xvalslong,yvalslong)
  
ax.set_xlim([-sep/2, sep/2])
ax.set_ylim([-sep/2, sep/2])
ax.set_zlim([-sep/2, sep/2])
ax.set_ylabel("y")
ax.set_xlabel("x")
ax.set_zlabel("z")
plt.show()

