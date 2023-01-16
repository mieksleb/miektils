#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:19:40 2022

2D plectoneme

@author: michaelselby
"""
import numpy as np
from tools import sn,cn,dn,E,F,am,K,E
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings
import tikzplotlib

warnings.filterwarnings("ignore")
# np.seterr(all="ignore")

r0 = 1

m_spec = 1/fsolve(lambda m: K(m)-2*E(np.pi/2,m),0.8261)[0]

def theta(s, theta0, m, A, f):
    lam = (A/f)**0.5
    phi = (theta0-np.pi)/2
    k = m**0.5
    theta = np.pi + 2*am(s/(lam*k)+F(phi,m),m)
    return theta

def theta_decrease(s, theta0, m, A, f):
    lam = (A/f)**0.5
    theta = theta0 - 2*am(s/(lam*m**0.5),m)
    return theta

def xfunc(s,theta0,thetaL,m,A,f):
    lam = (A/f)**0.5
    k = m**0.5
    theta = np.pi + 2*am(s/(lam*k) + F((theta0-np.pi)/2,m),m)
    xval = (2/m-1)*s - 2*lam/k*( E((theta-np.pi)/2,m)- E((theta0-np.pi)/2,m))
    return xval

def yfunc(s,theta0,thetaL,m,A,f):
    lam = (A/f)**0.5
    k = m**0.5
    u = s/(lam*k) + F((theta0-np.pi)/2,m)
    u0 = F((theta0-np.pi)/2,m)
    yval = 2 * lam/k * ( dn(u, m) - dn(u0, m))
    return yval

def x_max(theta0, m, A, f):
    lam = (A/f)**0.5
    val = (2-m) * ( F(np.pi/4,m)-F((np.pi-theta0)/2,m) ) - 2 * ( E(np.pi/4,m)-E((np.pi-theta0)/2,m) )
    val *= - lam/m**0.5
    return val

def xfunc_end(theta0,thetaL,m,A,f):
    lam = (A/f)**0.5
    k = m**0.5
    xval = (2-m)*(F((thetaL-np.pi)/2,m) - F((theta0-np.pi)/2,m)) - 2 * ( E((thetaL-np.pi)/2,m)- E((theta0-np.pi)/2,m))
    xval *=  lam / k
    return xval


class Elastica:
    def __init__(self, theta0, thetaL, L, origin, m, A, f, npoints):
        self.theta0 = theta0
        self.thetaL = thetaL
        self.L = L
        self.origin=origin
        self.m = m
        self.A = A
        self.f = f
        self.npoints = npoints
        self.lam = (A/f)**0.5
        
    def generate_elastica(self):
        k = self.m**0.5
        self.s = np.linspace(0,self.L,self.npoints)
        ds = self.s[1] - self.s[0]
        phi0 = np.full(self.npoints,(self.theta0-np.pi)/2)
        self.theta = np.pi + 2*am(self.s/(self.lam*k)+ F(phi0,self.m), self.m)
        # discrete version
        # self.x = np.cumsum(np.cos(self.theta)) * ds 
        # self.y = np.cumsum(np.sin(self.theta)) * ds 
        self.x = xfunc(self.s,self.theta0,self.thetaL,self.m,self.A,self.f)
        self.y = yfunc(self.s,self.theta0,self.thetaL,self.m,self.A,self.f)
        self.x += self.origin[0]
        self.y += self.origin[1]
        
    def disc_energy(self):
        self.ene = energy_discrete(self.theta, self.A, self.s, self.f)
        return self.ene
    
class Circular_Arc:
    def __init__(self,theta0, thetaL, L, origin, A, npoints):
        self.theta0 = theta0
        self.thetaL = thetaL
        self.L = L
        self.origin=origin
        self.A = A
        self.npoints = npoints
        
    def generate_arc(self):
        self.s = np.linspace(0, self.L, self.npoints)
        self.theta = self.theta0 + (self.thetaL-self.theta0)*self.s/self.L

        # 
        self.x = self.L * np.sin(self.theta) / (self.thetaL- self.theta0)
        self.y =  - self.L * np.cos(self.theta) / (self.thetaL- self.theta0)
        self.x -= self.L * np.sin(self.theta0) / (self.thetaL- self.theta0)
        self.y -= self.L * np.cos(self.theta0) / (self.thetaL- self.theta0)
        
    def disc_energy(self):
        self.ene = energy_discrete(self.theta, self.A, self.s, 0)
        return self.ene

        
        
    
def x_close(theta0,thetaL,m):
    phi0 = (theta0 - np.pi)/2
    phiL = (thetaL - np.pi)/2
    return (2-m)*(F(phiL,m)-F(phi0,m))-2*(E(phiL,m)-E(phi0,m))
  
def length( theta0, thetaL, m, A, f ):
    lam = (A/f)**0.5
    phi0 = (theta0 - np.pi)/2
    phiL = (thetaL - np.pi)/2
    l = F(phiL,m) - F(phi0,m)
    l *= lam*m**0.5
    return l

def energy_discrete(theta, A, s, f ):
    '''
    Calculates the discrete approximation to the energy of the loop.
    theta is a numpy array of theta values.
    A is the bending stiffness at each arclength value s

    '''
    ds = s[1]-s[0]
    n = len(s)
    diff = np.array([(theta[k+1]-theta[k])**2 for k in range(0,n-1)])
    bend = np.sum(diff)
    bend *= A/(2*ds)
    pull = np.sum(np.cos(theta))
    pull *= -f*ds
    return bend + pull
    


class Plectoneme:
    def __init__(self,L,A,f,nbraids,npoints,damage=False,damage_type="kinked",thetadam = 30*np.pi/180,A1=1,L1=1, m1=1):
        self.L = L
        self.A = A
        self.f = f
        self.nbraids = nbraids
        self.npoints = npoints
        self.damage = damage
        self.damage_type = damage_type
        self.thetadam = thetadam
        self.A1 = A1
        self.L1 = L1
        self.lam = (A/f)**0.5
        self.lam1 = (A1/f)**0.5
        self.m1 = m1
              
        
    def generate_plectoneme(self):
        mvals = np.linspace(0.01,2,50)
        lowest_ene = 10000
        for m_loop in mvals:
  
            if self.damage:
                if self.damage_type=="kinked":
                    thetaloopL = np.pi - self.thetadam/2
                    thetaloop0 = fsolve(lambda theta: x_close(theta,thetaloopL,m_loop),0.3)[0]
                    if x_close(thetaloop0,thetaloopL,m_loop) > 0.01:
                        continue
                    Lloop = length( thetaloop0, thetaloopL, m_loop, self.A, self.f)
                elif self.damage_type=="floppy":
                    thetaflop = np.pi - 2*am(self.L1/(2*self.lam1*self.m1**0.5),self.m1) 
                    thetaloop0 = fsolve(lambda theta: xfunc_end(thetaflop,np.pi,self.m1,self.A1,self.f) + xfunc_end(theta,thetaflop,m_loop,self.A,self.f),0.8)[0]
                    # print(xfunc_end(thetaflop,np.pi,m1,self.A1,self.f) + xfunc_end(thetaloop0,thetaflop,m_loop,self.A,self.f))
                    if thetaloop0 > np.pi/2:
                        continue
                    thetaloopL = thetaflop
                    # print(thetaloopL, "\n")
                    Lloop = length(thetaloop0,thetaloopL,m_loop,self.A,self.f)
                    # print(Lloop)
                    
                else:
                    print("Unrecognised Damage Type. Kinked or Floppy currently supported.")
                    print("If Kinked, enter a bending angle for thetadam.")
                    print("If Floppy, enter a value for A1 and a damage length L1.")
            else:
                thetaloopL = np.pi
                thetaloop0 = fsolve(lambda theta0: x_close(theta0,thetaloopL,m_loop),0.8)[0]
                Lloop = length(thetaloop0,thetaloopL,m_loop,self.A,self.f)
                
            if x_max(thetaloop0,m_loop,self.A,self.f) < r0:
                continue
            
            # once we have constructed the loop, the braid and tail follow using same process for damage/undamaged
    
            thetabraidL = np.pi - thetaloop0
            thetabraid0 = np.pi - thetabraidL
            Lbraid =  r0 * (thetabraidL - thetabraid0) / (1 - np.sin (thetabraid0))
            thetatailL = thetabraidL
            Ltail = self.L / 2 - self.nbraids * Lbraid - Lloop
            
            ene = 0
            # penalise configurations which don't have room for tail
            if Ltail < 0:
                ene += 10000
            
            if self.damage and self.damage_type=="floppy":
                Ltail -= self.L1/2
                ene += (1-2/self.m1)*self.f*self.L1/2 + 4*(self.A1*self.f)**0.5*(E((np.pi-thetaflop)/2,self.m1))
                
             # for long enough plectonemes m is approximately one, with energy difference minimal
            if Ltail/self.lam > 2: 
                # print("scree")
                m_tail = 1
            else:
                m_tail = fsolve(lambda m: Ltail - length(0,thetatailL,m,self.A,self.f), 0.999)[0]
            # calculate the energy of the plectoneme
            ene += (2/m_loop-1)*self.f*Lloop
            ene += self.nbraids*self.A * (thetabraidL - thetabraid0)**2 / (2 * Lbraid)
            # ene += 2*np.pi**2*self.nbraids**2 * self.A * np.sin(thetabraid0)**2 /Lbraid             
            ene += (1-2/m_tail)*self.f*Ltail + 4*(self.A*self.f)**0.5*(E(np.pi/2,m_tail)-E((thetatailL)/2,m_tail))
            if self.damage and self.damage_type=="kinked":
                print(thetaloop0)
                print(Lbraid)
                print(Lloop)
                print(Ltail , "\n")
            if ene < lowest_ene:

                self.Lloop = Lloop
                self.m_loop = m_loop
                self.thetabraid0 = thetabraid0
                self.thetabraidL = thetabraidL
                self.thetaloop0 = thetaloop0
                self.thetaloopL = thetaloopL
                self.thetatailL = thetatailL
                self.Ltail = Ltail
                self.Lbraid = Lbraid
                self.m_tail = m_tail
                self.exact_energy = ene
                lowest_ene = ene
                if self.damage and self.damage_type=="floppy":
                    self.thetaflop = thetaflop

        
        if lowest_ene >= 10000:
            print("Plectoneme not found")
        
        if self.damage and self.damage_type=="floppy":
            self.dam = Elastica(np.pi, self.thetaflop, self.L1/2, np.array([0,0]), self.m1, self.A, self.f, self.npoints)
            self.dam.generate_elastica()

        self.loop = Elastica(self.thetaloop0, self.thetaloopL , self.Lloop, np.array([0,0]), self.m_loop, self.A, self.f, self.npoints)
        self.loop.generate_elastica()
            
        
        self.braidlist = []
        for i in range(0,self.nbraids):
            braid = Circular_Arc(self.thetabraid0, self.thetabraidL, self.Lbraid, np.array([0,0]), self.A, self.npoints)
            braid.generate_arc()
            self.braidlist.append(braid)
            self.braidlist[i].x -= self.braidlist[i].x[-1]
            self.braidlist[i].y -= self.braidlist[i].y[-1]
            if i > 0:
                self.braidlist[i].y += self.braidlist[i-1].y[0]

            
        self.tail = Elastica(np.pi + self.thetatailL, 2*np.pi, self.Ltail, np.array([0,0]), self.m_tail, self.A, self.f, self.npoints)
        self.tail.generate_elastica()

        self.tail.x -= self.tail.x[0]
        self.tail.y -= self.tail.y[0]
        self.tail.y += self.braidlist[-1].y[0]
        
        y_min = self.tail.y[-1]
        for curve in [self.tail, self.loop] + self.braidlist:
            curve.y -= y_min
        if self.damage and self.damage_type=="floppy":
            self.dam.y -= self.dam.y[-1]
            self.dam.y += self.loop.y[-1]

        
        
        
    def calculate_energy(self,discrete=False):
        """
        Calculates half the energy of a plectoneme

        """
            
        if discrete == False:
            loop_ene = (2/self.loop.m-1)*self.f*self.loop.L 
            braid_ene = 2*np.pi**2*self.nbraids**2 * self.A * np.sin(self.thetabraid0)**2 / self.Lbraid
            tail_ene =  (1 - 2/self.tail.m) * self.f * self.tail.L + 4*(self.A*self.f)**0.5*(E(np.pi/2,self.tail.m)-E((self.thetatailL)/2,self.tail.m))/self.m_tail**0.5
        else:
            loop_ene = self.loop.disc_energy()
            tail_ene = self.tail.disc_energy()   
            braid_ene = self.nbraids * self.braidlist[0].disc_energy()
            
        ene = loop_ene + braid_ene + tail_ene
        
        if self.damage and self.damage_type=="floppy":
            if discrete:
                dam_ene = self.dam.disc_energy()
            else:
                dam_ene = (1-2/self.m1)*self.f*self.L1/2 + 4*(self.A1*self.f)**0.5*(E((np.pi-self.thetaflop)/2,self.m1))
                
            ene += dam_ene
        
        self.energy = ene
        return ene
        
    
    
    def print_values(self):
        if self.damage:
            if self.damage_type=="floppy":
                plec_type = "Floppy"
            else:
                plec_type = "Kinked"
        else:
            plec_type = "Undamaged"
                
            
        print( plec_type + " Plectoneme of length " +str(self.L)+ " with "+str(self.nbraids)+" braids.")
        print("Energy:", self.energy)
        
        
    def plot(self,tikz=False):
        plt.axis('equal')
        for i in range(0,self.nbraids):
            plt.plot(self.braidlist[i].x,self.braidlist[i].y, color='blue')
            plt.plot(-self.braidlist[i].x,self.braidlist[i].y, color='blue')
            
        plt.plot(self.tail.x,self.tail.y, color='blue', label='tail')
        plt.plot(self.loop.x,self.loop.y, color='blue',label='loop')
        plt.plot(-self.tail.x,self.tail.y, color='blue', label='tail')
        plt.plot(-self.loop.x,self.loop.y, color='blue',label='loop')
        
        if self.damage and self.damage_type=="floppy":
            plt.plot(self.dam.x,self.dam.y, color='red', label='dam')
            plt.plot(-self.dam.x,self.dam.y, color='red', label='dam')
            
        
        if tikz:
            tikzplotlib.save("plectoneme.tex")
        plt.show()
            
        
kT = 4.114  # kT in pN nm at 298K
A = 45  # bend persistence length in units of nm kbT        
        
homo = Plectoneme(L=200,A=A,f=2,nbraids=3,npoints=100)
homo.generate_plectoneme()
homo.calculate_energy(discrete=True)
homo.calculate_energy()
homo.plot(tikz=False)

print(" ")
        
kinked = Plectoneme(L=200,A=50,f=2,nbraids=3,npoints=100,damage=True,damage_type="kinked",thetadam = 30*np.pi/180)
kinked.generate_plectoneme()
kinked.calculate_energy(discrete=True)
kinked.plot()

lowest_energy = 10000
for m1 in np.linspace(0.8,2,50):
    floppy1 = Plectoneme(L=200,A=50,f=2,nbraids=3,npoints=100,damage=True,damage_type="floppy",A1 = 31,L1=3, m1=m1)
    floppy1.generate_plectoneme()
    floppy1.calculate_energy()
    energy = floppy1.energy
    if energy < lowest_energy:
        floppy = floppy1
floppy.plot()

homo.print_values()
floppy.print_values()
kinked.print_values()

