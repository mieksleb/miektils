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
L = 0.3


gump = scipy.special.ellipeinc(1, 1.52)

x = L*f**0.5/(2*(a1)**0.5)
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
    return 4*(a1*f)**0.5*(1/k)*(2*scipy.special.ellipeinc(am(u,k),k**2)-u)+8*(a*f)**0.5*(1-sn(u,k))+2*L*f

def Energy_Recip(alpha):
    return 8*(a1*f)**0.5*(scipy.special.ellipeinc(am(x,alpha),alpha**2)-x+0.5*x*alpha**2)+8*(a*f)**0.5*(1-alpha*sn(x,alpha))+2*L*f
    

def Energy(k):
    u=x/k
    return 4*(a1*f)**0.5*(1/k)*(2*k**2*sn(u,k)*cn(u,k)/dn(u,k)-x/k)+8*(a*f)**0.5*(1-sn(u,k))+2*L*f


# this function's zero is the minimum of the elastic energy, i.e it is the first derivative of the energy
def function(k):
    u = x/k
    return k*cn(u,k)/dn(u,k)-r**0.5

# another formulation of the energy, this only has dependence on k in the sn and sn^2 terms
def Energy2(k):
    u = x/k
    return 8*(a*f)**0.5*(1-sn(u,k)) + 8*(a1*f)**0.5*sn(u,k)+2*L*f

# linear approximation to the energy in the case that x is small
def Energy_Linear(k):
    return 8*(a*f)**0.5 + 2*L*f*(1-(1/r))

def Energy_Quadratic(k):
    return 8*(a*f)**0.5 + 2*L*f*(1-(1/r))-(1/6)*(L**3*a*f**2/(a1**2))*(1-(1/r))

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



fvals=np.linspace(0,10)