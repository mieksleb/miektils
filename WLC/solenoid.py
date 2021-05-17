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

A = 7
C = 20
bp = 100
lpbp = 3 #length per base pair in nm
L = bp*lpbp
kappa = 2
tau = 1
h = 10.5    # helical repeat in bp
Lk0 = bp/h  # Lk0 is resting linking number
wr0 = tau*((1+(kappa/tau)**2)**0.5-1) #mean field writhe
p = abs(tau)/(tau**2+kappa**2)
r = kappa/(tau**2+kappa**2)

def r(s): 
    theta = s/(r**2+p**2)**0.5
    return np.array([r*np.cos(theta),r*np.sin(theta),p*theta])

def tangent(s):
    theta = s/(r**2+p**2)**0.5
    return np.array([-r*np.sin(theta),r*np.cos(theta),p])/(r**2+p**2)**0.5

def integrand(s,t):
    chonk = np.dot(r(s)-r(t),r(s)-r(t))
    return np.dot(np.cross(t(s),t(t)),r(s)-r(t))/chonk**1.5


Wrnum = scipy.integrate.dblquad(integrand, 0, L, 0, L)
Wrnum = Wrnum//(2*np.pi*L)