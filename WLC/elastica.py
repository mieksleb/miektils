#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 10:29:47 2020

@author: michaelselby
"""
import scipy
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate



L=10     # length of the elastica
L1=0.5    # length of the floppy segnment
f=1     # force applied
A=1     # bending modulous of elastica
A1=0.5  # bending modulous of floppy segment elastica

lam1=(A1/f)**0.5
lam=(A/f)**0.5
r=lam1/lam
npoints=100
smin=-L/2
smax=L/2
s2min=-L1/2
s2max=L1/2


def sech(x):
    return (np.cosh(x))**(-1)


def h(lam,lam1,L1):
    return 2*lam1*(1-sech(L1/(lam1)))-2*lam*(1-sech(L1/(lam)))


s_flop = np.linspace(smin,smax,npoints)      # arclength range for inhomogeneity
s = np.linspace(s2min,s2max,npoints)   # arclength range for posiitve homogenous region

x_flop = [None]*npoints                      # x values for floppy segment
y_flop = [None]*npoints                      # y values for floppy segment
for i in range(0,npoints):
    y_flop[i] = 2*lam1*((sech((s_flop[i]/lam1)))-1)
    x_flop[i] = s_flop[i]-2*lam*np.tanh(s_flop[i]/lam1)
    
x = [None]*npoints
y = [None]*npoints
for i in range(0,npoints):
    y[i] = 2*lam1*(sech((s[i]/lam))-1)+h(lam,lam1,L1)
    x[i] = s[i]-2*lam*np.tanh(s[i]/lam)
    
del_list=[]
end_x = x_flop[0]
end_y = y_flop[0]
for i in range(0,npoints):
    if (-end_x<x[i]<end_x) and (y[i]>end_y):
        del_list.append(i)
        x[i]=0
        y[i]=0

        
print(x_flop[0],y_flop[0])
print(x_flop[-1],y_flop[-1])

xx = [None]*del_list[0]         # x values for hom loop with floppy region removed
yy = [None]*del_list[0]         # y values for hom loop with floppy region removed
for i in range(0,del_list[0]):
    xx[i] = x[i]
    yy[i] = y[i]
    
xx_flop = [None]*(npoints-del_list[-1]-1)
yy_flop = [None]*(npoints-del_list[-1]-1)
for i in range(0,npoints-del_list[-1]-1):
    xx_flop[i] = x[del_list[-1]+i+1]
    yy_flop[i] = y[del_list[-1]+i+1]
    
    
plt.plot(xx,yy, color='blue')
plt.plot(-xx,yy,color='blue')
plt.plot(x_flop[0],y_flop[0],color='red',marker='x')
plt.plot(-x_flop[0],y_flop[0],color='red',marker='x')
plt.plot(xx_flop,yy_flop,color='orange')
plt.show()
  
'''
We now calculate the energy of the loop using the integral formula for the elastica hamiltonian

'''




    