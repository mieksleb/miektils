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
from scipy.integrate import quad
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

kT = 300*1.38064852*10**-23

A = 150 #bend persistence length in units of bp kbT
C = 300
f = 1
bp = 200
L = bp  #length in bp
L_braid = 40
L_tail = 10
w = 0.5
npoints = 200
wpoints = 2
Lk = -1
p = 3.4
Lk0 = bp/p
lam = (A/f)**0.5
u0 = 0.9
phi0 = 0
psi0 = 0
L_loop = 34.44948
 


'''
Define all ellitpic functions

'''
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
def K(m):
    return scipy.special.ellipk(m)
def E_comp(m):
    return scipy.special.ellipe(m)

def Pi(n,m):
    func = lambda theta : 1 / ( (1-n*np.sin(theta)**2) * (1 - m * np.sin(theta)**2)**0.5 )
    val = scipy.integrate.quad(func, 0, np.pi/2)[0]
    return val

def Pi_inc(n,u,m):
    func = lambda theta : 1 / ( (1-n*np.sin(theta)**2) * (1 - m * np.sin(theta)**2)**0.5 )
    val = scipy.integrate.quad(func, 0, u)[0]
    return val



# define the rotation matrix R1
def rot1(theta):
    return np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
# defines the rotation matrix R3
def rot3(theta):
    return np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
# defines the rotation R3(psi)*R1(theta)*R3(phi)
def euler_rot(psi,theta,phi):
    return np.matmul(np.matmul(rot3(psi),rot1(theta)),rot3(phi))

svals = np.linspace(0,L/2,npoints)






def u0_function(s,lam):
    return s-2*lam*np.tanh(s/lam)


# s_int is the value of s where the homoclinic elasticae self intersects
s_int = fsolve(lambda s: u0_function(s,lam),-L_loop)[0]

# u0 is then the cosine of the angle at the point of intersection
u0 = s_int**2/(2*lam**2) -1
theta0 = np.arccos(u0)*180/np.pi
sintheta0 = (1-u0**2)**(0.5)


class elastica:
    def __init__(self, A, C, f, L, a, b, c, m, npoints):
        self.A = A
        self.C = C
        self.f = f
        self.L = L
        self.a = a
        self.b = b
        self.c = c
        self.u0 = u0
        self.m = m

    def lam(self):
        return (self.A / self.f)**0.5
        
    def svals(self):
        return np.linspace(0,self.L/2,npoints)
    
    def p1(self):
        return (2 * self.A * self.f * (self.a + 1) * (1 + self.b) * (1 + self.c))**0.5
    def p2(self):
        return (2 * self.A * self.f * (self.a - 1) * (1 - self.b) * (1 - self.c))**0.5
   
    def p_psi(self):
        return (self.p1() + self.p2())/2
    
    def p_phi(self):
        return (self.p1() - self.p2())/2
    
    def gamma(self):
        return self.a + self.b + self.c - (self.p_psi())**2/(2 * self.A * self.f)
    
    def energy(self):
        pull = (self.gamma()/2 - self.a) * self.f * self.L
        bend = 2 * ( 2 * self.A * self.f * (self.a - self.c))**0.5 * E(am(self.L / (2 * self.lam() * (self.m **0.5)), self.m), self.m)
        return pull + bend
    
    def costheta(self,s):
        return (self.b - self.c) * ( sn( s * ((self.a - self.c)/2)**0.5 / self.lam(), self.m ))**2 + self.c
    
    def curvature(self,s):
        return 2*(self.gamma()-self.costheta(self,s))/(self.lam())**2

    def theta(self,s):
        u = self.costheta(s)

        val = np.arccos(u)
        # if s < 0: 
        #     val += np.pi
        return val
    
    def phi(self,s):
        u = s * ((self.a - self.c)/2)**0.5/self.lam()
        factor = 1 /(2 * self.A * self.f * (self.a - self.c))**0.5
        if np.isclose(self.c,-1) == True:
            term1 = 0
        else:
            term1 = (self.p_phi()+self.p_psi()/(self.c+1)) * Pi_inc((self.c-self.b)/(self.c+1),u,self.m)
        if np.isclose(self.p_psi()-self.p_phi(),0):
            term2 = 0
        else:
            term2 = (self.p_psi()-self.p_phi())/(self.c-1) * Pi_inc((self.b-self.c)/(self.c-1),u,self.m)
        val = factor * (term1 - term2)
        return val
    
    def psi(self,s):
        return self.p_psi()*(1/self.C-1/self.A)*s + self.phi(s)
        
    
    def thetavals(self):    
        thetavals = []
        for s in self.svals():
            thetavals.append(self.theta(s))
        return thetavals
    
    def phivals(self):    
        phivals = []
        for s in self.svals():
            phivals.append(self.phi(s))
        return phivals
    
    def psivals(self):    
        psivals = []
        for s in self.svals():
            psivals.append(self.psi(s))
        return psivals
    
    def printvals(self):
        print("a = "+ str(self.a))
        print("b = "+ str(self.b))
        print("c = "+ str(self.c))
        print("p_psi = "+ str(self.p_psi()))
        print("p_phi = "+ str(self.p_phi()))
        print("gamma = "+ str(self.gamma()))
        print("Energy = "+ str(self.energy())+"\n")
        
class complex_elastica:
    def __init__(self, A, C, f, L, rea, ima, b, m, constA, p_psi1, npoints):
        self.A = A
        self.C = C
        self.f = f
        self.L = L
        self.rea = rea
        self.ima = ima
        self.b = b
        self.m = m
        self.constA = constA
        self.p_psi1 = p_psi1

    def lam(self):
        return (self.A / self.f)**0.5
    
    # def constA(self):
    #     return (self.rea-self.b)**2+self.ima**2
        
    def svals(self):
        return np.linspace(-self.L/2,self.L/2,npoints)
    
    # def p1(self):
    #     return (2 * self.A * self.f * ((1+self.rea)**2+self.ima**2) * (1 + self.b))**0.5
    # def p2(self):
    #     return (2 * self.A * self.f * ((1-self.rea)**2+self.ima**2) * (1 - self.b))**0.5
   
    def p_psi(self):
        return self.p_psi1
    
    def p_phi(self):
        return -self.p_psi1
    
    def gamma(self):
        return np.real(2*self.rea + self.b - (self.p_psi())**2/(2 * self.A * self.f))
    
    def q(self):
        return (self.constA-2)/(self.constA+2)
    
    # def energy(self):
    #     pull = (self.gamma()/2 - self.a) * self.f * self.L
    #     bend = 2*( 2 * self.A * self.f * (self.a - self.c))**0.5 * E(am(self.L / (2 * self.lam() * (self.m **0.5)), self.m), self.m)
    #     return pull + bend
    
    def costheta(self,s):
        u = s * (2*self.constA)**0.5/self.lam()
        return self.b + self.constA * np.tan(am(u,self.m)/2)**2
    
    def curvature(self,s):
        return 2*(self.gamma()-self.costheta(self,s))/(self.lam())**2

    def theta(self,s):
        u = self.costheta(s)
        val = np.arccos(u)
        # if u < 0: 
        #       val += np.pi
        return val
    
    def phi(self,s):
        u = s * (2*self.constA**0.5/self.lam())
        # print(am(u,self.m))
        factor1 = self.p_phi()/self.A
        term1 = s/(2+self.constA)

        qp = 1 - self.q()**2
        factor2 =   2 * self.constA * ( 1 - 1/ (1 - self.q()**2) ) / (self.constA+2)

        expr = (self.m + self.q()**2*(1-self.m))**0.5
        term2 = - self.q()**2 * Pi_inc( 1/qp , am(u,self.m), self.m)/qp
        term3 = self.q() * np.arctanh( sn(u,self.m) * expr / ( dn(u,self.m) * (1 - self.q()**2)**0.5)) / ( (1-self.q()**2)**0.5 * expr)
        val = factor1 * (term1 + factor2 * (term2 + term3))
        return val
    
    def psi(self,s):
        return self.p_psi()*(1/self.C-1/self.A)*s + 2*self.phi(s)
    
    def thetavals(self):    
        thetavals = []
        for s in self.svals():
            thetavals.append(self.theta(s))
        return thetavals
    
    def phivals(self):    
        phivals = []
        for s in self.svals():
            phivals.append(self.phi(s))
        return phivals
    
    def psivals(self):    
        psivals = []
        for s in self.svals():
            psivals.append(self.psi(s))
        return psivals
    
    def printvals(self):
        print("Re(a) = "+ str(self.rea))
        print("Im(a) = "+ str(self.ima))
        print("b = "+ str(self.b))
        print("m = "+ str(self.m))
        print("Constant A = "+ str(self.constA))
        print("p_psi = "+ str(self.p_psi()))
        print("p_phi = "+ str(self.p_phi()))
        print("gamma = "+ str(self.gamma()))
        print("q = "+ str(self.q()))
        # print("Energy = "+ str(self.energy())+"\n")

def euler2cartesian(elastica,triad0,origin):
    svals = elastica.svals()
    thetavals = elastica.thetavals()
    phivals = elastica.phivals()
    psivals = elastica.psivals()
    tangent_vals = []
    normal_vals = []
    binormal_vals = []
    tangent = triad0[0]
    normal = triad0[1]
    binormal = triad0[2]
    tangent_vals.append(tangent)
    normal_vals.append(normal)
    binormal_vals.append(binormal)
        
    for i in range(0,len(thetavals)):
        tangent = np.dot(euler_rot(psivals[i],thetavals[i],phivals[i]),tangent0)
        normal = np.dot(euler_rot(psivals[i],thetavals[i],phivals[i]),normal0)
        binormal = np.dot(euler_rot(psivals[i],thetavals[i],phivals[i]),binormal0)
        tangent_vals.append(tangent)
        normal_vals.append(normal)
        binormal_vals.append(binormal)
       
    ds = abs(svals[0]-svals[-1])/(len(svals)-1)
        
    rvals = []  
    r = origin
    rvals.append(r)
    for i in range(1,len(svals)):
        r = r + tangent_vals[i]*ds   # position is integral of the tangent
        rvals.append(r)
            
    return rvals



'''
Solving the Restraints

The type of constraints depend on what type of solution we are searching for

For closed loops, f is a lagrenge multiplier which must be determined, whereas for MT loops, it is an input

'''



##                      Loopus Baybee                        ##
###############################################################



def function(a,b,c):
    m = (b-c)/(a-c)
    return a*K(m)-(a-c)*E_comp(m)

def curv_function(a_new,a,b,c):
    return a_new*(1+b**2) + (a_new**2-1)**0.5*(1-b**2) - 2*b - (a-1)/(1+b)

def loop_function(b,L,m,lam):
    return L/(2*lam*m**0.5)-(2/(b+1))**0.5*K(m)

def loop_function2(b,m):
    return E_comp(m)/K(m)-1+m/(1+b)

# def m_loop(a,b):
#     return (b+1)/(a+1)




def theta_loop(s,b,m,lam):
    u = s*((b+1)/2)**0.5/(lam*m**0.5)
    val = np.arccos((b+1)*(sn(u,m))**2-1)
    if u < 0: 
        val *= -1
    return val 


def phi_loop(s,a,b,m):
    u = (s*((b+1)/2)**0.5)/(lam*m**0.5)
    return ((1-b)*(a-1)/(2*(a+1)))**0.5*Pi_inc(-(b+1)/2,u,m)





def theta_homoclinic(s,lam,m):
    u = s/lam
    return np.pi+2*am(u,m)

    





'''
Generating the values

'''

origin = np.array([0,0,0])
tangent0 = np.array([0,0,1])
normal0 = np.array([0,1,0])
binormal0 = np.cross(tangent0,normal0)
triad0 = [tangent0,normal0,binormal0]

r = L/lam

def poundland(r,m):
    K1 = K(m)
    alpha = r**2/(4*m*K1**2)
    return r/np.pi * (1 + 1/(m*(alpha-1)))**0.5 * (Pi(1/(1-alpha),m)/K(m) - 0.5 * (1-1/alpha))-1

print(Pi(1/(1-r**2/(4*0.85*K(0.85**2))),0.85))
print("chook")

m_loop = fsolve(lambda m: poundland(r,m)-poundland(r,0.95),0.96)[0]
print(m_loop)

chumvals = []
mvals = np.linspace(0,1,50)
for m in mvals:
    chumvals.append(poundland(r,m))
    
# plt.plot(mvals,chumvals)
# plt.show()

L_loop = 100
m_loop = 0.8
# Loop 
c_loop = -1
b_loop = 0.5
a_loop = 1.8

p = (A*f)**0.5*(a_loop + 1) * (c_loop + 1)




#generate the loop
loop = elastica(A, C, f, L_loop, a_loop, b_loop, c_loop, m_loop, npoints)
# loop = elastica(A, C, f, L, 1, 1, -1, 1, 100)
loop.printvals()
ds = svals[1]-svals[0]
svals = loop.svals()
rvals = euler2cartesian(loop, triad0, origin)

tangent0 = np.array([0,0,-1])
normal0 = np.array([0,-1,0])
binormal0 = np.cross(tangent0,normal0)
triad0 = [tangent0,normal0,binormal0]
rvals2 = euler2cartesian(loop, triad0, origin)






'''
Plotting

'''


fig = plt.figure()
ax = fig.gca(projection='3d')

xvals = []
yvals = []
zvals = []
xvals2 = []
yvals2 = []
zvals2 = []


for j in range(0,len(rvals)):
        xvals.append(rvals[j][0])  
        yvals.append(rvals[j][1])      
        zvals.append(rvals[j][2])
        xvals2.append(rvals2[j][0])  
        yvals2.append(-rvals2[j][1])      
        zvals2.append(rvals2[j][2])
      

        
        
# xsep  = max(xvals)-min(xvals)
# ysep  = max(yvals)-min(yvals)
# zsep  = max(zvals)-min(zvals)

# sep = max([xsep,ysep,zsep])

xvalslong = []
yvalslong = []
zvalslong = []


        
ax.scatter(zvals,xvals,yvals, c='red',s=0.1)
ax.scatter(zvals2,xvals2,yvals2, c='blue',s=0.1)

# ax.scatter(zvals3,xvals3,yvals3, c='green',s=0.1)
# ax.scatter(zvals4,xvals4,yvals4, c='purple',s=0.1)
# ax.scatter(zvals5,xvals5,yvals5, c='orange',s=0.1)
# ax.scatter(0,0,0, c='black',s=1)
# ax.plot3D(zvalslong,xvalslong,yvalslong)
  
# ax.set_xlim([-sep, sep])
# ax.set_ylim([-sep, sep])
# ax.set_zlim([-sep, sep])
ax.set_ylabel("y")
ax.set_xlabel("x")
ax.set_zlabel("z")
plt.show()

