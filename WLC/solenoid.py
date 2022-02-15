#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 13:24:38 2021


Solenoidz

@author: michaelselby
"""

import scipy
import numpy as np
from scipy import integrate
import numba
from numba import njit
from numba import jit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

T = 300      # temperature in Kelvin
kT = 1.38064852*10**(-23)*T  # kT 
beta = 1/kT
A = 7*10**(-9)*kT  # bending modulus 
C = 20*10**(-9)*kT
f = 50*10**(-12)
lam = (A/f)**0.5
bp = 100
lpbp = 3*10**(-10) #length per base pair in m
L = bp*lpbp
p = 20*lpbp
radius = 10*lpbp

kappa = radius/(radius**2+p**2)
tau = p/(radius**2+p**2)
cosdelta = 1/(1+(kappa/tau)**2)**0.5
mu = cosdelta/lam**2
h = 10.5    # helical repeat in bp
Lk0 = bp/h  # Lk0 is resting linking number
sigma = 0.05
deltaLk = sigma*Lk0
Lk = Lk0+deltaLk
wr0 = tau*((1+(kappa/tau)**2)**0.5-1)*L/(2*np.pi) #mean field writhe
p = abs(tau)/(tau**2+kappa**2)
radius = kappa/(tau**2+kappa**2)
npoints = 100


def r(s): 
    theta = s/(radius**2+p**2)**0.5
    return np.array([radius*np.cos(theta),radius*np.sin(theta),p*theta])

def tangent(s):
    theta = s/(radius**2+p**2)**0.5
    return np.array([-radius*np.sin(theta),radius*np.cos(theta),p])/(radius**2+p**2)**0.5

def integrand(s,t):
    chonk = np.dot(r(s)-r(t),r(s)-r(t))
    return np.dot(np.cross(tangent(s),tangent(t)),r(s)-r(t))/chonk**1.5


svals = np.linspace(0,L,npoints)
bpivals = np.linspace(0,bp,npoints)
    
    
rvals = []
for s in svals:
    rvals.append(r(svals))
    
xx = [vec[0] for vec in rvals]
yy = [vec[1] for vec in rvals]
zz = [vec[2] for vec in rvals]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(xx,yy,zz, c='red',s=0.1)
# ax.plot3D(xx,yy,zz)
plt.show()


def get_spline(veclist,smin,smax, k = 3, s = 0, per = False):
   
    import scipy.interpolate
 
    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in veclist]
    yy = [vec[1] for vec in veclist]
    zz = [vec[2] for vec in veclist]
    
    vals = np.linspace(smin,smax,len(xx))    

    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(vals, xx, k = 3, per = per)
    spline_yy = scipy.interpolate.splrep(vals, yy, k = 3, per = per)
    spline_zz = scipy.interpolate.splrep(vals, zz, k = 3 , per = per)
    
    
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

@njit(fastmath=True)
def norm(vec):
    return vec / np.sqrt(np.dot(vec,vec))

# @njit(fastmath=True)
def get_solenoid_writhe(L,pitch,radius, npoints = 100, circular = False, integral_type = "simple"):
    

    import scipy.interpolate
    
    svals = np.linspace(0,L,npoints)
    
    
    rvals = []
    for s in svals:
        rvals.append(r(s))
        
        
        
    spline = get_spline(rvals,0,L)
    

    s1xx, s1yy, s1zz = spline[0]

    

    xx = scipy.interpolate.splev(svals, s1xx)
    yy = scipy.interpolate.splev(svals, s1yy)
    zz = scipy.interpolate.splev(svals, s1zz)
    

    dmxx = scipy.interpolate.splev(svals, s1xx, 1)
    dmyy = scipy.interpolate.splev(svals, s1yy, 1)
    dmzz = scipy.interpolate.splev(svals, s1zz, 1)

    tt = list(range(len(svals)))
    for ii in range(len(svals)):
        tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])
  
    
    # get the normal vector via n(s) = dt(s)/ds = d^2/ds^2[r(s)]    
    ddmxx = scipy.interpolate.splev(svals, s1xx, 2)
    ddmyy = scipy.interpolate.splev(svals, s1yy, 2)
    ddmzz = scipy.interpolate.splev(svals, s1zz, 2)
    nn = list(range(len(svals)))
    for ii in range(len(svals)):
        nn[ii] = np.array([ddmxx[ii], ddmyy[ii], ddmzz[ii]])
 
    
    ds = L/len(svals)
    # do the integration w.r.t. s
 
    
    # triple vector product using numba fastmast
    @njit(fastmath=True)
    def multidet(u,v,w):
        d=np.empty(3)
        d=\
        u[0]*(v[1]*w[2]-v[2]*w[1])+\
        u[1]*(v[2]*w[0]-v[0]*w[2])+\
        u[2]*(v[0]*w[1]-v[1]*w[0])  # 14 operations / det
        return d

    
    
    # now for the integration which we use numba fastmath capabilities for speed
    @njit(fastmath=True)
    def discrete_dbl_int(tt, xx, yy, zz, dx, dy, ss, circular = circular):
        if circular:
            srange = range(len(ss)-1)
        else:
            srange = range(len(ss))

        writhe_integral = 0
        for ii in srange:
            for jj in srange:
                # skip ii=jj and use symmetry in {ii, jj}
                if ii > jj:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = multidet((tt[ii]), tt[jj], diff_frac)
                    writhe_integral += triple_scalar_product
            
        writhe_integral *= dx * dy / (2 * np.pi) # multiplied by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
        return writhe_integral
    
    
    writhe = discrete_dbl_int(tt, xx, yy, zz, ds, ds, svals, circular = circular)
    
    return writhe
    
#writhe = get_solenoid_writhe(L,p,radius, npoints = 100, circular = False, integral_type = "simple")


'''
We now calculate the Free Energy, and tangent autocorrelation function in both real and Fourier space 
'''


omega = 2*np.pi**2*C*(deltaLk-wr0)/(A*L)

def tan_corr_fourier(q):
    alpha = q**2+cosdelta/lam**2+tau**2-omega*tau
    beta = ((2*tau-omega)*q)**2
    return (2*alpha-kappa**2)/(alpha*(alpha-kappa**2)-beta)

qvals = np.linspace(-10,10,npoints)

tan_corr_fourier_vals = []
for q in qvals:
    tan_corr_fourier_vals.append(tan_corr_fourier(q))
    
plt.plot(qvals,tan_corr_fourier_vals)
plt.show()


bee = 2*mu+2*tau**2-2*omega*tau-kappa**2-(2*tau-omega)**2
cee = (mu+tau**2-omega*tau)*(mu+tau**2-omega*tau-kappa**2)

a1 = mu+tau**2-omega*tau
b1 = (2*tau-omega)**2
c1 = kappa**2

gamma = (-a1/2+(c1+b1)/4+0.5*(a1**2-a1*c1)**0.5)**0.5
eta = (a1/2-(c1+b1)/4+0.5*(a1**2-a1*c1)**0.5)**0.5




z = (-bee/2+((bee/2)**2-cee)**0.5)**0.5
# gamma = z.real
# eta = z.imag

Free_energy = eta*L

def tan_corr(s):
    return 1-eta*L**2/(beta*A)+eta*L**2*np.exp(-eta*abs(s))/(beta*A)*(np.cos(gamma*s)-(eta/gamma)*np.sin(gamma*s))

tan_corr_vals = []
for s in svals:
    tan_corr_vals.append(tan_corr(s))
    
plt.plot(bpivals,tan_corr_vals)
plt.show()

E_bend  = A*L*kappa**2
E_twist = 2*np.pi*C*deltaLk**2/L
E_ext = -f*L

Etot = (E_bend+E_twist+E_ext)*beta

print(eta*L**2/(beta*A))
