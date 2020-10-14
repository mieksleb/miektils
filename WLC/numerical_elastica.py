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
from scipy.optimize import curve_fit


# parameters
# a1 is the bending modulus of the floppy segment
# L is the length of the floppy segment
# a is the bending modulus of the rest of the curve
# f is the applied force

a1 = 0.4
a = 1
f = 1
L1 = 0.3
L = 10
theta0=np.pi/10
npoints=100



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

# define the incomplete elliptic integral of the second kind
def E(theta,k):
    return scipy.special.ellipeinc(theta,k**2)

def F(theta,k):
    return scipy.special.ellipkinc(theta,k**2)






def Energy_floppy(k,a,a1,L1,r,f):
    x = L1*f**0.5/(2*(a1)**0.5)
    u=x/k
    return 4*(a1*f)**0.5*(1/k)*(2*E(am(u,k),k)-u)+8*(a*f)**0.5*(1-sn(u,k))+2*L1*f

# this function's zero is the minimum of the elastic energy, i.e it is realted to  the first derivative of the energy
def function_floppy(k,a,a1,L1,r,f):
    x = L1*f**0.5/(2*(a1)**0.5)
    u = x/k
    return k*cn(u,k)/dn(u,k)-r**0.5


# linear approximation to the energy in the case that x is small
def Energy_Linear_floppy(a,a1,L1,r,f):
    r=a1/a
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/r))


def Energy_Quadratic_floppy(a,a1,L1,r,f):
    r=a1/a
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/r))-(1/6)*(L1**3*a*f**2/(a1**2))*(1-(1/r))


def Energy_Cubic_floppy(a,a1,L1,r,f):
    r=a1/a
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/r))-(1/6)*(L1**3*a*f**2/(a1**2))*(1-(1/r))+(1/72)*(L1**5*f**3)*(a**2/a1**4)*(1-(1/r))

def Energy_Homoclinic_floppy(a,a1,L1,r,f):
    return 8*(a*f)**0.5*(1-np.tanh(x))+8*(a1*f)**0.5*np.tanh(x)


def Energy_bent(a,a1,L1,r,f,theta0):
    return 8*(a*f)*(1-np.sin(theta0/4))



def Energy_bent_floppy(k,a,a1,L1,r,f,theta0):
    x = L1*f**0.5/(2*(a1)**0.5)
    u=x/k+F(theta0/4,k)
    return 4*(a1*f)**0.5*(1/k)*(2*E(am(u,k),k)-2*E(theta0/4,k)-x/k)+8*(a*f)**0.5*(1-sn(u,k))+2*L1*f


# this function's zero is the minimum of the elastic energy, i.e it is realted to  the first derivative of the energy
def function_bent_floppy(k,a,a1,L1,r,f,theta0):
    x = L1*f**0.5/(2*(a1)**0.5)
    u = x/k
    return k*cn(u,k)/dn(u,k)-r**0.5


# linear approximation to the energy in the case that x is small
def Energy_Linear_bent_floppy(k,a,a1,L1,r,f,theta0):
    r=a1/a
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/r)) - 2*(a1*f)**0.5*theta0

    

# this function tests whether two functions, fun1 and fun2, are equivalent on some interval with values
# array, to within some tolerance tol. Output is True if equivalent and False if not.
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




k_num = minimize(lambda k: Energy_floppy(k,a,a1,L1,r,f),r**0.5).x[0]
E_num = minimize(lambda k: Energy_floppy(k,a,a1,L1,r,f),r**0.5).fun
k_exact = fsolve(lambda k: function_floppy(k,a,a1,L1,r,f),r**0.5)[0]
Energy_Exact = Energy_floppy(k_exact,a,a1,L1,r,f)
k_lin = r**0.5
E_lin = Energy_Linear_floppy(a,a1,L1,r,f)
E_quad = Energy_Quadratic_floppy(a,a1,L1,r,f)
E_cub = Energy_Cubic_floppy(a,a1,L1,r,f)
E_homo = Energy_Homoclinic_floppy(a,a1,L1,r,f)



def Energy_bongoboi(c,k,a,a1,L1,r,f):
    return np.abs(8*(a*f)**0.5+2*L1*f*(1-(1/r))*1/(1+(a*L1**2*f/(c*a1**2)))-Energy_floppy(k,a,a1,L1,r,f))
    



a1vals = np.linspace(0.1,1,25)

c_vals=[]
for a1val in a1vals:
    a1 = a1val
    fun = lambda k: function_floppy(k,a,a1,L1,r,f)
    r=a1/a
    k_exact = fsolve(fun,r**0.5)[0]
    Energy_Exact = Energy_floppy(k_exact,a,a1,L1,r,f)
    fun2 = lambda c: Energy_bongoboi(c,k_exact,a,a1,L1,r,f)
    c_choi = minimize(fun2,15).x[0]
    c_vals.append(c_choi) 


def func(x,m,n):
    return 12+m/(x**n)

opt = curve_fit(func, a1vals, c_vals)

mopt=opt[0][0]
nopt=opt[0][1]

a1vals = a1vals[:-1]
c_vals = c_vals[:-1]


c_theory_vals =[] 
for a1val in a1vals:
    c_theory_vals.append(func(a1val, mopt, nopt))

plt.plot(a1vals,c_vals)
plt.plot(a1vals,c_theory_vals)
plt.show()






k_bent_flop_num = minimize(lambda k: Energy_bent_floppy(k,a,a1,L1,r,f,theta0),r**0.5).x[0]
E_bent_flop_num = minimize(lambda k: Energy_bent_floppy(k,a,a1,L1,r,f,theta0),r**0.5).fun
E_bent_flop = Energy_bent_floppy(k_exact,a,a1,L1,r,f,theta0)
E_bent_flop_lin = Energy_Linear_bent_floppy(k_exact,a,a1,L1,r,f,theta0)


print("Elliptic Modulus (Exact)", k_exact)
print("Elliptic Modulus (Numerical)", k_num)
print("Elliptic Modulus (Linear)", k_lin)
print("\n")
print("Energy (Exact)", Energy_Exact)
print("Energy (Numerical)", E_num)
print("Energy (Linear)", E_lin)
print("Energy (Quadratic)", E_quad)
print("Energy (Homoclinic)", E_homo)
print("\n")

print("Elliptic Modulus Bent (Numerical)", k_bent_flop_num)
print("\n")
print("Energy Bent Flop (Exact)", E_bent_flop)
print("Energy Bent Flop (Numerical)", E_bent_flop_num)
print("Energy Bent Flop (Linear)", E_bent_flop_lin)





# TEP data
tenbp = [3.333,6.389,10,0]
fourbp = [3.611,4.722,6.111,5.278]
twobp = [1.667,2.222,1.778,3.333]
onebp = [2.778,2.222,1.389,1.667]
force = [1,2,3,4]



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

def x(k,L,lam,s):
    return s*((2/k**2)-1)-(2*lam/k)*scipy.special.ellipeinc(am(s/(k*lam),k),k**2)

def y(k,L,lam,s):
    return (2*lam/k)*dn(s/(k*lam),k)


sright = np.linspace(L1/2,L/2,npoints)
sleft = np.linspace(-L/2,-L1/2,npoints)
s1 = np.linspace(-L1/2,L1/2,npoints)


plt.plot(x(0.8,L1,lam1,s1),y(0.8,L1,lam1,s1)+(y(1,L1,lam,L1/2)-y(0.8,L1,lam1,L1/2)),color='red')
plt.plot(x(1,L,lam,sleft),y(1,L,lam,sleft),color='blue')
plt.plot(x(1,L,lam,sright),y(1,L,lam,sright),color='blue')
plt.show()

