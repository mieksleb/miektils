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
import tikzplotlib


# parameters
# KT is kBT in nmpN
# a1 is the bending modulus of the floppy segment with units 
# L is the length of the floppy segment in nm
# a is the bending modulus of the rest of the curve
# f is the applied force in pN


KT = 4.28
a = 45*KT
a1 = 0.7*a
f = 2
L1 = 3.3
L = 100
theta0=30*np.pi/180
npoints=100



x = L1*f**0.5/(2*(a1)**0.5)
print(x)

r = a1/a
print(r)
rootr = r**0.5
print(rootr)
 
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


def Energy_floppy_simple(k,a,a1,L1,r,f):
    return 8*(a*f)**0.5 + 2*L1*f*(1-(1/k)**2)

# this function's zero is the minimum of the elastic energy, i.e it is realted to  the first derivative of the energy
def function_floppy(k,a,a1,L1,r,f):
    x = L1*f**0.5/(2*(a1)**0.5)
    u = x/k
    return k*cn(u,k)/dn(u,k)-r**0.5




def function_iter(a,a1,L1,r,f):
    x = L1*f**0.5/(2*(a1)**0.5)
    r=a1/a
    k=r**0.5
    nsteps=50
    for i in range(0,nsteps):
        k = r**0.5/np.cos(x*(1-r)**0.5/k)
    return k



def functionnewton_raphson(a,a1,L1,r,f):
    x = L1*f**0.5/(2*(a1)**0.5)
    r=a1/a
    alpha=1
    nsteps=50
    for i in range(0,nsteps):
        alpha = alpha - (alpha+np.cos(alpha*x*((1/r)-1)**0.5))/(1+np.sin(alpha*x*((1/r)-1)**0.5))
    return -r**0.5/alpha
    

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
    return 8*(a*f)**0.5*(1-np.sin(theta0/4))


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



E_bent = Energy_bent(a, a1, L1, r, f, theta0)
E_bent_lin = 8*(a*f)**0.5*(1-theta0/4)


k_num = minimize(lambda k: Energy_floppy(k,a,a1,L1,r,f),r**0.5).x[0]
E_num = minimize(lambda k: Energy_floppy(k,a,a1,L1,r,f),r**0.5).fun
k_exact = fsolve(lambda k: function_floppy(k,a,a1,L1,r,f),r**0.5)[0]
k_iter = function_iter(a, a1, L1, r, f)
k_newton = functionnewton_raphson(a, a1, L1, r, f)
Energy_Exact = Energy_floppy(k_exact,a,a1,L1,r,f)
Energy_Exact_simple = Energy_floppy_simple(k_exact,a,a1,L1,r,f)
k_lin = r**0.5
k_quad = r**0.5+(1-r)/(2*r**0.5)*x**2
k_quart = r**0.5+(1-r)*x**2/(2*r**0.5)+(r-1)*(r+7)*x**4/(24*r**1.5)
k_sext = r**0.5+(1-r)*x**2/(2*r**0.5)+(r-1)*(r+7)*x**4/(24*r**1.5)-(r-1)*(0.17)*x**6/(r**2.5)
E_lin = Energy_Linear_floppy(a,a1,L1,r,f)
E_quad = Energy_Quadratic_floppy(a,a1,L1,r,f)
E_quart = 8*(a*f)**0.5 + 2*L1*f*(1-(1/r))-(1/6)*(L1**3*a*f**2/(a1**2))*(1-(1/r))-(1-(1/r))*2*L1*f*x**4*(2*(r-2)/(15*r**2))
E_cub = Energy_Cubic_floppy(a,a1,L1,r,f)
E_homo = Energy_Homoclinic_floppy(a,a1,L1,r,f)

E_flop_diff = 8*a**0.5*f**0.5-Energy_Exact
E_flop_diff_lin = Energy_Exact-E_lin
E_diff_lin  = 2*L1*f*(1-1/r)
E_flop_diff_quad = Energy_Exact-E_quad


k_bent_flop_num = minimize(lambda k: Energy_bent_floppy(k,a,a1,L1,r,f,theta0),r**0.5).x[0]
E_bent_flop_num = minimize(lambda k: Energy_bent_floppy(k,a,a1,L1,r,f,theta0),r**0.5).fun
E_bent_flop = Energy_bent_floppy(k_exact,a,a1,L1,r,f,theta0)
E_bent_flop_lin = Energy_Linear_bent_floppy(k_exact,a,a1,L1,r,f,theta0)

print("-----------------------------")
print("Bent Chain")
print("\n")
print("Energy (Exact)", E_bent)
print("Energy Bent (Linear)", E_bent_lin)
print("-----------------------------")



print("-----------------------------")
print("Floppy Chain")
print("\n")
print("Elliptic Modulus (Exact)", k_exact)
print("Elliptic Modulus (Numerical)", k_num)
print("Elliptic Modulus (Linear)", k_lin)
print("Elliptic Modulus (Quadratic)", k_quad)
print("Elliptic Modulus (Quartic)", k_quart)
print("\n")
print("Energy (Exact)", Energy_Exact)
print("Energy (Numerical)", E_num)
print("Energy (Linear)", E_lin)
print("Energy (Quadratic)", E_quad)
print("Energy (Homoclinic)", E_homo)
print("-----------------------------")
print("\n")


print("-----------------------------")
print("Bent and Floppy Chain")
print("\n")
print("Elliptic Modulus (Numerical)", k_bent_flop_num)
print("\n")
print("Energy (Exact)", E_bent_flop)
print("Energy (Numerical)", E_bent_flop_num)
print("Energy (Linear)", E_bent_flop_lin)
print("-----------------------------")





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



'''
plot of the bent sequence
'''

L_eff = 2*lam*np.log(np.tan(((np.pi+theta0/2)/4)))
sright2 = np.linspace(L_eff/2,L/2,npoints)
sleft2 = np.linspace(-L/2,-L_eff/2,npoints)
xright2 = x(1,L,lam,sright2[-1])
plt.plot(x(1,L,lam,sleft2)-x(1,L,lam,sleft2[-1]),y(1,L,lam,sleft2),color='blue')
plt.plot(x(1,L,lam,sright2)-x(1,L,lam,sright2[0]),y(1,L,lam,sright2),color='blue')
plt.xlim([-xright2, xright2])
plt.ylim([-10, 15])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('bent')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
# tikzplotlib.save("bent.tex")




'''
plot of the floppy sequence
'''

sright = np.linspace(L1/2,L/2,npoints)
sleft = np.linspace(-L/2,-L1/2,npoints)
s1 = np.linspace(-L1/2,L1/2,npoints)

xright = x(1,L,lam,sright[-1])
plt.plot(x(k_exact,L1,lam1,s1),y(k_exact,L1,lam1,s1)+(y(1,L1,lam,L1/2)-y(k_exact,L1,lam1,L1/2)),color='blue')
plt.plot(x(1,L,lam,sleft)+(x(1,L1,lam,L1/2)-x(k_exact,L1,lam1,L1/2)),y(1,L,lam,sleft),color='blue')
plt.plot(x(1,L,lam,sright)-(x(1,L1,lam,L1/2)-x(k_exact,L1,lam1,L1/2)),y(1,L,lam,sright),color='blue')
plt.xlim([-xright, xright])
# plt.xlim([-20, 20])
plt.ylim([-10, 15])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('floppy')
plt.xlabel('x')
plt.ylabel('y')
# tikzplotlib.save("flop.tex")


'''
plot of the bent and floppy sequence
'''
L_eff2 = 2*lam1*np.log(np.tan(((np.pi+theta0/2)/4)))
sright = np.linspace(L1/2,L/2,npoints)
sleft = np.linspace(-L/2,-L1/2,npoints)
s2 = np.linspace(-L1/2,-L_eff2/2,npoints)
s3 = np.linspace(L_eff2/2,L1/2,npoints)

xright = x(1,L,lam,sright[-1])
plt.plot(x(k_exact,L1,lam1,s2)-x(k_exact,L,lam1,s2[-1]),y(k_exact,L1,lam1,s2)+(y(1,L1,lam,L1/2)-y(k_exact,L1,lam1,L1/2)),color='red')
plt.plot(x(k_exact,L1,lam1,s3)-x(k_exact,L,lam1,s3[0]),y(k_exact,L1,lam1,s3)+(y(1,L1,lam,L1/2)-y(k_exact,L1,lam1,L1/2)),color='red')
plt.plot(x(1,L,lam,sleft)+(x(1,L1,lam,L1/2)-x(k_exact,L1,lam1,L1/2))-x(k_exact,L,lam1,s2[-1]),y(1,L,lam,sleft),color='blue')
plt.plot(x(1,L,lam,sright)-(x(1,L1,lam,L1/2)-x(k_exact,L1,lam1,L1/2))-x(k_exact,L,lam1,s3[0]),y(1,L,lam,sright),color='blue')
plt.xlim([-xright, xright])
plt.ylim([-10, 15])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('bent and floppy')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
# tikzplotlib.save("bent_flop.tex")




def E_boi(a,a1,L1,r,f,k):
    x = L1*f**0.5/(2*(a1)**0.5)
    print(r**0.5/k)
    theta = np.arcsin(r**0.5/k)
    print(theta)
    return (E(theta,k)-E(np.pi/2,k))-(r*(1-r))**0.5/((k**2-r)**0.5)

smol = E_boi(a,a1,L1,r,f,k_exact)
print("smol",smol)





x = L1*f**0.5/(2*(a1)**0.5)
u = x/k_exact
k = k_exact
m=k**2

honk = F(np.arcsin(r**0.5/k),k)-F(np.pi/2,k)+x/k
print("honk",honk)


kprime2 = 1-k_exact**2
theta = np.arcsin(r**0.5/k_exact)
chob = -r**0.5/((k**2)*((1-r)**0.5)*(1-(r/k**2))**0.5)+(E(theta,k))/(k*kprime2)-(F(theta,k)/k)-(r/(1-r))**0.5*(1-(r/k**2))**0.5/kprime2
print(chob)
chib = x/(k**2)+E(np.pi/2,k)/(k*kprime2)-F(np.pi/2,k)/k
print(chib)

print(chob-chib)

