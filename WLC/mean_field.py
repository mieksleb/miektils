#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:07:02 2020

@author: michaelselby
"""
      
import sys
import base
import numpy as np
import tools
import readers

import matplotlib.pyplot as plt
from numba import njit
from numba import jit
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt
import scipy
import cmath
import numba
import scipy.interpolate



npoints=100

conf = "miekbin/circle_300-Lk7.dat"
top = "miekbin/circle_300-Lk7.top"

circular=True

l = readers.LorenzoReader(conf, top)
s = l.get_system()

strand1 = s._strands[0]
strand2 = s._strands[1]

bp = len(strand1._nucleotides)

spline1 = tools.get_bb_spline(strand1)
spline2 = tools.get_bb_spline(strand2, reverse = True)   

s1xx, s1yy, s1zz = spline1[0]
s2xx, s2yy, s2zz = spline2[0]

smin = spline1[1][0]
smax = spline1[1][-1]

# bpi is the base pair index parameter that is common to both splines
bpi = np.linspace(smin, smax, npoints)

# find the midpoint between the input splines, as a function of base pair index
m1xx = (scipy.interpolate.splev(bpi, s1xx) + scipy.interpolate.splev(bpi, s2xx)) / 2
m1yy = (scipy.interpolate.splev(bpi, s1yy) + scipy.interpolate.splev(bpi, s2yy)) / 2
m1zz = (scipy.interpolate.splev(bpi, s1zz) + scipy.interpolate.splev(bpi, s2zz)) / 2


# contour_len[ii] is contour length along the midpoint curve of point ii
delta_s = [np.sqrt((m1xx[ii+1]-m1xx[ii])**2+(m1yy[ii+1]-m1yy[ii])**2+(m1zz[ii+1]-m1zz[ii])**2) for ii in range(len(bpi)-1)]
contour_len = np.cumsum(delta_s)
contour_len = np.insert(contour_len, 0, 0)

# ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
ss = np.linspace(contour_len[0], contour_len[-1], npoints)

# get the splines as a function of contour length
msxx = scipy.interpolate.splrep(contour_len, m1xx, k = 3, s = 0, per = circular)
msyy = scipy.interpolate.splrep(contour_len, m1yy, k = 3, s = 0, per = circular)
mszz = scipy.interpolate.splrep(contour_len, m1zz, k = 3, s = 0, per = circular)


xx = scipy.interpolate.splev(ss, msxx)
yy = scipy.interpolate.splev(ss, msyy)
zz = scipy.interpolate.splev(ss, mszz)



# find the tangent of the midpoint spline. 
# the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
# get discrete array of normalised tangent vectors; __call__(xxx, 1) returns the first derivative
# the tangent vector is a unit vector
dmxx = scipy.interpolate.splev(ss, msxx, 1)
dmyy = scipy.interpolate.splev(ss, msyy, 1)
dmzz = scipy.interpolate.splev(ss, mszz, 1)




# we also need the 'normal' vector u(s) which points between the base pairs. (or between the spline fits through the backbones in this case)
# n.b. these uxx, uyy, uzz are not normalised
uxx_bpi = scipy.interpolate.splev(bpi, s2xx) - scipy.interpolate.splev(bpi, s1xx)
uyy_bpi = scipy.interpolate.splev(bpi, s2yy) - scipy.interpolate.splev(bpi, s1yy)
uzz_bpi = scipy.interpolate.splev(bpi, s2zz) - scipy.interpolate.splev(bpi, s1zz)

# get the normal vector spline as a function of contour length
suxx = scipy.interpolate.splrep(contour_len, uxx_bpi, k = 3, s = 0, per = circular)
suyy = scipy.interpolate.splrep(contour_len, uyy_bpi, k = 3, s = 0, per = circular)
suzz = scipy.interpolate.splrep(contour_len, uzz_bpi, k = 3, s = 0, per = circular)
    


dogs = [[suxx[i], suyy[i], suzz[i]] for i in range(len(xx))]

tt = list(range(len(ss)))
for ii in range(len(ss)):
    tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])
  

ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
# do the integration w.r.t. s
 



