#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:27:42 2020

@author: michaelselby
"""
import scipy
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt


# parameters
# a1 is the bending modulus of the floppy segment
# L is the length of the floppy segment
# a is the bending modulus of the rest of the curve
# f is the applied force

a1 = 0.65
a = 1
f = 1
L1 = 0.3
L = 10



x = L1*f**0.5/(2*(a1)**0.5)
r = a1/a


# define the 4 Jacobi elliptic functions: sn, cn, dn and am

def sn(u,k):
    return scipy.special.ellipj(u,k**2)[0]

def cn(u,k):
    return scipy.special.ellipj(u,k**2)[1]

def dn(u,k):
    return scipy.special.ellipj(u,k**2)[2]

def am(u,k):
    return scipy.special.ellipj(u,k**2)[3]

def Energy_Full(k):
    u=x/k
    return 4*(a1*f)**0.5*(1/k)*(2*scipy.special.ellipeinc(am(u,k),k**2)-u)+8*(a*f)**0.5*(1-sn(u,k))+2*L1*f

def Energy_Recip(alpha):
    return 8*(a1*f)**0.5*(scipy.special.ellipeinc(am(x,alpha),alpha**2)-x+0.5*x*alpha**2)+8*(a*f)**0.5*(1-alpha*sn(x,alpha))+2*L1*f
    

def Energy(k):
    u=x/k
    return 4*(a1*f)**0.5*(1/k)*(2*k**2*sn(u,k)*cn(u,k)/dn(u,k)-x/k)+8*(a*f)**0.5*(1-sn(u,k))+2*L1*f


# this function's zero is the minimum of the elastic energy, i.e it is the first derivative of the energy
def function(k):
    u = x/k
    return k*cn(u,k)/dn(u,k)-r**0.5

# another formulation of the energy, this only has dependence on k in the sn and sn^2 terms
def Energy2(k):
    u = x/k
    return 8*(a*f)**0.5*(1-sn(u,k)) + 8*(a1*f)**0.5*sn(u,k)+2*L1*f

# linear approximation to the energy in the case that x is small
def Energy_Linear(k):
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/r))

def Energy_Quadratic(k):
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/r))-(1/6)*(L1**3*a*f**2/(a1**2))*(1-(1/r))

def Energy_Homoclinic(k):
    return 8*(a*f)**0.5*(1-np.tanh(x))+8*(a1*f)**0.5*np.tanh(x)
    

# this function tests whether two functions, fun1 and fun2, are equivalent on some interval with values
# array, to within some tolerance tol. Output is True is equivalent and False if not.
def are_equivalent(fun1,fun2,array,tol):
    per = True
    for val in array:
        if np.abs(fun1(val)-fun2(val)) < tol:
            per = True
            continue
        else:
            per = False
            break
        print(per)
    return per




k_num = minimize(Energy_Full,r**0.5).x[0]
E_num = minimize(Energy_Full,r**0.5).fun

k_recip = 1/(minimize(Energy_Recip,r**0.5).x[0])
E_recip = minimize(Energy_Recip,r**0.5).fun[0]


k_exact = fsolve(function,r**0.5)[0]
Energy_exact = Energy(k_exact)
Energy_full = Energy_Full(k_exact)
Energy_dog = Energy2(k_exact)

k_lin = r**0.5
E_lin = Energy_Linear(k_exact)
E_quad = Energy_Quadratic(k_exact)
E_homo = Energy_Homoclinic(k_exact)

test = are_equivalent(Energy_Full, Energy, np.linspace(0.01,1), 0.1)
k_vals = np.linspace(0.3,1,100)
E1vals = []
E2vals = []
for val in k_vals:
    E1vals.append(Energy_Full(val))
    E2vals.append(Energy(val))

plt.plot(k_vals,E1vals)
plt.plot(k_vals,E2vals)
plt.show()


tonk = Energy_Recip(1/k_exact)


print("Elliptic Modulus (Exact)", k_exact)
print("Elliptic Modulus (Numerical)", k_num)
print("Elliptic Modulus (Linear)", k_lin)

print("\n")

print("Energy (full)", Energy_full)
print("Energy (Exact)", Energy_exact)
print("Energy (Numerical)", E_num)
print("Energy (Linear)", E_lin)
print("Energy (Quadratic)", E_quad)
print("Energy (Homoclinic)", E_homo)





'''
Now for the plots
'''


lam1=(a1/f)**0.5
lam=(a/f)**0.5
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
s = np.linspace(s2min,s2max,npoints)         # arclength range for posiitve homogenous region

x_flop = [None]*npoints                      # x values for floppy segment
y_flop = [None]*npoints                      # y values for floppy segment
for i in range(0,npoints):
    y_flop[i] = 2*lam1*((sech((s_flop[i]/lam1)))-1)
    x_flop[i] = s_flop[i]-2*lam*np.tanh(s_flop[i]/lam1)
    
x = [None]*npoints
y = [None]*npoints
for i in range(0,npoints):
    y[i] = 2*lam1*(sech((s[i]/lam))-1)+h(lam,lam1,L)
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
  







